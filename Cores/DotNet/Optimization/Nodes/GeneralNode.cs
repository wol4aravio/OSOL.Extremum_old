﻿using System;
using System.Collections.Generic;

namespace OSOL.Extremum.Cores.DotNet.Optimization.Nodes
{
    
    using Area = Dictionary<string, Tuple<double, double>>;

    public abstract class GeneralNode<TBase, TFuncType, TV> where TV: class, IOptimizable<TBase, TFuncType>
    {
        public int NodeId;
        
        public abstract void Initialize(Func<Dictionary<string, TFuncType>, TFuncType> f, Area area, State<TBase, TFuncType, TV> state);
        public abstract void Process(Func<Dictionary<string, TFuncType>, TFuncType> f, Area area, State<TBase, TFuncType, TV> state);

        public virtual int? GetCurrentCondition(Func<Dictionary<string, TFuncType>, TFuncType> f, Area area, State<TBase, TFuncType, TV> state) => null;
    }
}