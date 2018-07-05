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
        private static string currentSolutionEfficiencyName = "currentSolutionEfficiency"; // f_x

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
                state.SetParameter(currentSolutionEfficiencyName, x.GetPerformance(f));
            }

            public override void Process(Func<Dictionary<string, double>, double> f, Area area, State<RealVector, double, RealVector> state)
            {

            }
            
        }
          
        private class AdaptiveSearchNode : GeneralNode<RealVector, double, RealVector>
        {

            public AdaptiveSearchNode(int nodeId)
            {
                this.NodeId = nodeId;
            }

            public override void Initialize(Func<Dictionary<string, double>, double> f, Area area, State<RealVector, double, RealVector> state)
            {
                
            }

            public override void Process(Func<Dictionary<string, double>, double> f, Area area, State<RealVector, double, RealVector> state)
            {
                var x = state.GetParameter<RealVector>(currentSolutionName);
                var f_x = state.GetParameter<double>(currentSolutionEfficiencyName);

                var alpha = state.GetParameter<double>(expansionCoefficientName);
                var beta = state.GetParameter<double>(compressionCoefficientName);
                
                var M = state.GetParameter<int>(numberOfBadTrialsName);
                var N1 = state.GetParameter<int>(numberOfIterationsOfAdaptiveSearchName);
                
                var t_min = state.GetParameter<double>(minStepValueName);
                
                var t = 1.0;
                var k = 0;
                
                for (int j = 0; j < N1; ++j)
                {
                    var xi = gorn.GetContinuousUniformVector(area.ToDictionary(
                        kvp => kvp.Key, 
                        kvp => Tuple.Create(-1.0, 1.0)));
                    var xi_norm = Math.Sqrt(xi.Select(kvp => kvp.Value * kvp.Value).Sum());

                    var x_new = x.MoveBy(xi.ToDictionary(kvp => kvp.Key, kvp => t * kvp.Value / xi_norm));
                    x_new = x_new.Constrain(area);
                    
                    var f_x_new = x_new.GetPerformance(f);
                    if (f_x_new < f_x)
                    {
                        var z = x.MoveBy(((x_new - x) * alpha).Elements);
                        var f_z = z.GetPerformance(f);

                        x = x_new;

                        if (f_z < f_x_new)
                        {
                            x = z;
                            f_x = f_z;
                            t *= alpha;
                            continue;
                        }
                    }

                    k++;
                    if (k == M)
                    {
                        k = 0;
                        if (t < t_min)
                        {
                            return;
                        }
                        else
                        {
                            t *= beta;
                        }
                    }
                }
            }
            
        }

    }
}