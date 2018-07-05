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
        
        private static GoRN gorn = new GoRN();
        
    }
}