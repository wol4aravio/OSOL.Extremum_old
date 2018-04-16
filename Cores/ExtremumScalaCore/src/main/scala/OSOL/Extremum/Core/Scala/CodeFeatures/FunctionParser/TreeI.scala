package OSOL.Extremum.Core.Scala.CodeFeatures.FunctionParser

import OSOL.Extremum.Core.Scala.Vectors.VectorObject
import OSOL.Extremum.Core.Scala.Arithmetics.Interval

object TreeI {

  abstract class BinaryOpTree[Interval](left: Tree[Interval], right: Tree[Interval], op: (Interval, Interval) => Interval) extends Tree[Interval] {
    final override def calculate(v: VectorObject[Interval]): Interval = op(left.calculate(v), right.calculate(v))
  }

  abstract class UnaryOpTree[Interval](subTree: Tree[Interval], op: Interval => Interval) extends Tree[Interval] {
    final override def calculate(v: VectorObject[Interval]): Interval = op(subTree.calculate(v))
  }

  class AdditionTree(left: Tree[Interval], right: Tree[Interval]) extends BinaryOpTree[Interval](left, right, (x, y) => x + y)

  class SubtractionTree(left: Tree[Interval], right: Tree[Interval]) extends BinaryOpTree[Interval](left, right, (x, y) => x - y)

  class MultiplicationTree(left: Tree[Interval], right: Tree[Interval]) extends BinaryOpTree[Interval](left, right, (x, y) => x * y)

  class DivisionTree(left: Tree[Interval], right: Tree[Interval]) extends BinaryOpTree[Interval](left, right, (x, y) => x / y)

  class PowerTree(left: Tree[Interval], right: Tree[Interval]) extends BinaryOpTree[Interval](left, right, (x, y) => x ** y)

  class NegTree(subTree: Tree[Interval]) extends UnaryOpTree[Interval](subTree, x => -x)

  class SinTree(subTree: Tree[Interval]) extends UnaryOpTree[Interval](subTree, x => x.sin())

  class CosTree(subTree: Tree[Interval]) extends UnaryOpTree[Interval](subTree, x => x.cos())

  class ExpTree(subTree: Tree[Interval]) extends UnaryOpTree[Interval](subTree, x => x.exp())

  class AbsTree(subTree: Tree[Interval]) extends UnaryOpTree[Interval](subTree, x => x.abs())

  class LnTree(subTree: Tree[Interval]) extends UnaryOpTree[Interval](subTree, x => x.ln())

  class SqrtTree(subTree: Tree[Interval]) extends UnaryOpTree[Interval](subTree, x => x.sqrt())

  class ConstantTree(value: Interval) extends Tree[Interval] {
    final override def calculate(v: VectorObject[Interval]): Interval = value
  }

  class VariableTree(varName: String) extends Tree[Interval] {
    final override def calculate(v: VectorObject[Interval]): Interval = v(varName)
  }

}