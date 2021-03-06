package OSOL.Extremum.Apps.JVM

import OSOL.Extremum.Algorithms
import OSOL.Extremum.Cores.JVM.Optimization.Algorithm
import OSOL.Extremum.Cores.JVM.Optimization.RemoteFunctions.{IntervalRemoteFunction, RealRemoteFunction}
import OSOL.Extremum.Cores.JVM.Vectors.{IntervalVector, RealVector}
import OSOL.Extremum.Cores.JVM.Vectors.RealVector.Converters._
import org.rogach.scallop.ScallopConf
import spray.json._
import java.io._

import OSOL.Extremum.Cores.JVM.Arithmetics.Interval

object Runner extends App {

  class Conf(arguments: Seq[String]) extends ScallopConf(arguments) {
    val algorithm = opt[String](required = true)
    val task = opt[String](required = true)
    val port = opt[Int](required = true)
    val field = opt[String](default = Some("target"))
    val result = opt[String](required = true)
    val output = opt[String](default = Some("json"))
    val logStates = opt[String](default = None)
    val seed = opt[String](default = None)
    verify()
  }

  override def main(args: Array[String]): Unit = {
    val conf = new Conf(args)
    main(conf)
  }

  def runRealVectorAlgorithm(algorithm: Algorithm[RealVector, java.lang.Double, RealVector],
                             f: RealRemoteFunction,
                             area: Map[String, (java.lang.Double, java.lang.Double)], logStates: Option[String]): RealVector = {
    f.initialize()
    val result = algorithm.work(x => f(x), area, logStates)
    f.terminate()
    result
  }

  def runIntervalVectorAlgorithm(algorithm: Algorithm[IntervalVector, Interval, IntervalVector],
                                 f: IntervalRemoteFunction,
                                 area: Map[String, (java.lang.Double, java.lang.Double)], logStates: Option[String]): IntervalVector = {
    f.initialize()
    val result = algorithm.work(x => f(x), area, logStates)
    f.terminate()
    result
  }

  def saveRealVectorResult(result: RealVector, output: String, file: String): Unit = {
    output match {
      case "json" =>
        val out = new FileWriter(s"$file.json")
        out.write(result.convertToJson().prettyPrint)
        out.close()
      case "csv" =>
        val out = new FileWriter(s"$file.csv")
        out.write(result.elements.keys.mkString(",") + "\n")
        out.write(result.elements.values.mkString(",") + "\n")
        out.close()
    }
  }

  def saveIntervalVectorResult(result: IntervalVector, output: String, file: String): Unit = {
    output match {
      case "json" =>
        val out = new FileWriter(s"$file.json")
        out.write(result.toJson.prettyPrint)
        out.close()
      case "csv" =>
        val out = new FileWriter(s"$file.csv")
        out.write(result.elements.keys.mkString(",") + "\n")
        out.write(result.elements.values.mkString(",") + "\n")
        out.close()
    }
  }

  def main(conf: Conf): Unit = {
    val algConfig = scala.io.Source.fromFile(conf.algorithm()).getLines().mkString("\n").parseJson.asJsObject
    val Seq(language, algorithm) = algConfig
      .getFields("language", "algorithm")
      .map(j => j.asInstanceOf[JsString].value)

    val taskConfig = scala.io.Source.fromFile(conf.task()).getLines().mkString("\n").parseJson.asJsObject
    val area = taskConfig
      .getFields("area").head
      .asInstanceOf[JsArray].elements
      .map { j =>
        val Seq(JsString(name), JsNumber(min), JsNumber(max)) = j.asJsObject.getFields("name", "min", "max")
        (name, (new java.lang.Double(min.toDouble), new java.lang.Double(max.toDouble)))
      }.toMap[String, (java.lang.Double, java.lang.Double)]

    val logStates = conf.logStates.toOption

    (language, algorithm) match {
      case ("Scala", "RS") | ("Scala", "RandomSearch") => {
        val Seq(JsNumber(radius), JsNumber(maxTime)) = algConfig.getFields("radius", "maxTime")
        val algorithm =
          if (conf.seed.isEmpty) Algorithms.Scala.RandomSearch.createFixedStepRandomSearch(radius.toDouble, maxTime.toDouble)
          else {
            val seed = scala.io.Source.fromFile(conf.seed()).getLines().mkString("\n").parseJson.convertTo[RealVector]
            Algorithms.Scala.RandomSearch.createFixedStepRandomSearch(radius.toDouble, maxTime.toDouble, Some(seed))
          }

        val f = new RealRemoteFunction(conf.task(), conf.port(), conf.field())
        val result = runRealVectorAlgorithm(algorithm, f, area, logStates)

        saveRealVectorResult(result, conf.output(), conf.result())
      }
      case ("Java", "RS") | ("Java", "RandomSearch") => {
        val Seq(JsNumber(radius), JsNumber(maxTime)) = algConfig.getFields("radius", "maxTime")
        val algorithm = Algorithms.Java.RandomSearch.createFixedStepRandomSearch(radius.toDouble, maxTime.toDouble)

        val f = new RealRemoteFunction(conf.task(), conf.port(), conf.field())
        val result = runRealVectorAlgorithm(algorithm, f, area, logStates)

        saveRealVectorResult(result, conf.output(), conf.result())
      }
      case ("Scala", "IES") | ("Scala", "IntervalExplosionSearch") => {
        val Seq(JsNumber(maxBombs), JsNumber(rMaxRatio), JsNumber(maxTime)) = algConfig.getFields("maxBombs", "rMaxRatio", "maxTime")
        val algorithm =
          if (conf.seed.isEmpty) Algorithms.Scala.IntervalExplosionSearch.createIntervalExplosionSearch(maxBombs.toInt, rMaxRatio.toDouble, maxTime.toDouble)
          else {
            val seed: Seq[IntervalVector] = scala.io.Source.fromFile(conf.seed()).
              getLines()
              .mkString("\n")
              .parseJson
              .asInstanceOf[JsArray]
              .elements
              .map(_.convertTo[IntervalVector])
            Algorithms.Scala.IntervalExplosionSearch.createIntervalExplosionSearch(maxBombs.toInt, rMaxRatio.toDouble, maxTime.toDouble, Some(seed))
          }

        val f = new IntervalRemoteFunction(conf.task(), conf.port(), conf.field())
        val result = runIntervalVectorAlgorithm(algorithm, f, area, logStates)

        saveIntervalVectorResult(result, conf.output(), conf.result())
        saveRealVectorResult(result.toBasicForm().elements, conf.output(), conf.result() + "_real")

        val seedValues =
          new JsArray(Algorithms.Scala.IntervalExplosionSearch.extractSeed(algorithm).map(_.convertToJson()).toVector)
        val writer = new FileWriter(conf.result() + "_seed.json")
        writer.write(seedValues.prettyPrint)
        writer.close()

      }
      case _ => throw new Exception("Unsupported Algorithm")
    }
  }

}
