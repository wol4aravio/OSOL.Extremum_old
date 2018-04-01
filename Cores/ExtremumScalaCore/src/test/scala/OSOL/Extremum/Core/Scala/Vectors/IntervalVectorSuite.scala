package OSOL.Extremum.Core.Scala.Vectors

import OSOL.Extremum.Core.Scala.Arithmetics.Interval
import org.scalatest.FunSuite
import OSOL.Extremum.Core.Scala.Vectors.IntervalVector.Converters._
import OSOL.Extremum.Core.Scala.Vectors.Exceptions._
import OSOL.Extremum.Core.Scala.CodeFeatures.Pipe

class IntervalVectorSuite extends FunSuite {

  val v1: IntervalVector = Map("x" -> Interval(1.0), "y" -> Interval(2.0, 3.0), "z" -> Interval(3.0, 5.0))
  val v2: IntervalVector = Map("x" -> Interval(1.0), "z" -> Interval(-3.0, -2.0))

  test("Keys") {
    assert(v1.keys == Set("x", "y", "z"))
  }

  test("Value Extraction") {

    assert(v1("x") == Interval(1.0))
    assert(v1("y") == Interval(2.0, 3.0))
    assert(v1("z", Interval(0.0)) == Interval(3.0, 5.0))
    assert(v1("a", Interval(0.0)) == Interval(0.0))
    intercept[MissingKeyException] { v1("a") }

  }

  test("To String") {
    assert(v1.toString == "x -> [1.0; 1.0]\ny -> [2.0; 3.0]\nz -> [3.0; 5.0]")
  }

  test("Addition") {
    assert(v1 + v1 == v1 * 2.0)
    intercept[DifferentKeysException] { v1 + v2 }
  }

  test("Addition with Imputation") {
    assert(v1 ~+ v2 == (Map("x" -> Interval(2.0), "y" -> v1("y"), "z" -> Interval(0.0, 3.0)) |> IntervalVector.apply))
  }

  test("Subtraction") {
    assert(v1 - v1 == (Map("x" -> Interval(0.0), "y" -> Interval(-1.0, 1.0), "z" -> Interval(-2.0, 2.0)) |> IntervalVector.apply))
    intercept[DifferentKeysException] { v1 - v2 }
  }

  test("Subtraction with Imputation") {
    assert(v1 ~- v2 == (Map("x" -> Interval(0.0), "y" -> v1("y"), "z" -> Interval(5.0, 8.0)) |> IntervalVector.apply))
  }

  test("Multiplication") {
    assert(v1 * v1 == (Map("x" -> Interval(1.0), "y" -> Interval(4.0, 9.0), "z" -> Interval(9.0, 25.0)) |> IntervalVector.apply))
    intercept[DifferentKeysException] { v1 * v2 }

  }

  test("Multiplication with Imputation") {
    assert(v1 ~* v2 == (Map("x" -> Interval(1.0), "y" -> Interval(2.0, 3.0), "z" -> Interval(-15.0, -6.0)) |> IntervalVector.apply))
  }

  test("Multiply by coefficient") {
    assert(v1 + v1 == v1 * 2.0)
  }

  test("Move by") {
    assert(v1.moveBy("x" -> -1.0).moveBy(Seq("z" -> -3.0, "y" -> -2.0)) == (Map("x" -> Interval(0.0), "y" -> Interval(0.0, 1.0), "z" -> Interval(0.0, 2.0)) |> IntervalVector.apply))
  }

  test("Constraining") {
    assert(v1
      .constrain("x" -> (-1.0, 0.0))
      .constrain(Seq("y" -> (3.0, 10.0)))
      .constrain("z" -> (-5.0, 4.0)) == (Map("x" -> Interval(0.0), "y" -> Interval(3.0), "z" -> Interval(3.0, 4.0)) |> IntervalVector.apply))
  }

  test("Get Performance") {
    val f: Map[String, Interval] => Interval = v => v("x") - v("y") + v("z")
    assert(v1.getPerformance(f) == 1.0)
    assert((v1 * 2.0).getPerformance(f) == 2.0)
  }


}