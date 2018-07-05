using System;
using System.Collections.Generic;
using System.Linq;
using OSOL.Extremum.Cores.DotNet.Optimization;
using OSOL.Extremum.Cores.DotNet.Vectors;
using OSOL.Extremum.Cores.DotNet.Optimization.Nodes;
using OSOL.Extremum.Cores.DotNet.Random;
using OSOL.Extremum.Cores.DotNet.Random.Distributions;

namespace OSOL.Extremum.Algorithms.CSharp
{

    using Area = Dictionary<string, Tuple<double, double>>;

    public static class HybridRandomSearch
    {
        private static string numberOfGeneratedPointsName = "numberOfGeneratedPoints"; // R
        private static string numberOfAttemptsName = "numberOfAttempts"; // P
        private static string numberOfIterationsPerAttemptName = "numberOfIterationsPerAttempt"; // ITER
        private static string numberOfBadTrialsName = "numberOfBadTrials"; // M
        private static string numberOfIterationsOfAdaptiveSearchName = "numberOfIterationsOfAdaptiveSearch"; // N1

        private static string reductionCoefficientName = "reductionCoefficient"; // gamma
        private static string recoverCoefficientName = "recoverCoefficient"; // eta
        private static string expansionCoefficientName = "expansionCoefficient"; // alpha
        private static string compressionCoefficientName = "compressionCoefficient"; // beta

        private static string minStepValueName = "minStepValue"; // t_min
        private static string minDeltaForTargetFunctionName = "minDeltaForTargetFunction"; // eps_1
        private static string minDeltaForSolutionVectorName = "minDeltaForSolutionVector"; // eps_2

        private static string currentSolutionName = "currentSolution"; // x

        private static GoRN gorn = new GoRN();

        private class InitializationNode : GeneralNode<RealVector, double, RealVector>
        {

            public InitializationNode(int nodeId)
            {
                this.NodeId = nodeId;
            }

            public override void Initialize(Func<Dictionary<string, double>, double> f, Area area, State<RealVector, double, RealVector> state)
            {

                RealVector x = area.ToDictionary(kvp => kvp.Key, kvp => 0.5 * (kvp.Value.Item1 + kvp.Value.Item2));
                state.SetParameter(currentSolutionName, x);
            }

            public override void Process(Func<Dictionary<string, double>, double> f, Area area,
                State<RealVector, double, RealVector> state)
            {

            }
        }
    }
}