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
        private static string numberOfIterationsOfBestSampleName = "numberOfIterationsOfBestSample"; // N2s

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

        private static RealVector GenerateXi(Area area)
        {
            RealVector xi = gorn.GetContinuousUniformVector(area.ToDictionary(
                kvp => kvp.Key, 
                kvp => Tuple.Create(-1.0, 1.0)));

            return xi;
        }

        private static RealVector MoveByXi(RealVector x, RealVector xi, Area area, double t)
        {
            var xi_norm = Math.Sqrt(xi.Elements.Select(kvp => kvp.Value * kvp.Value).Sum());

            var x_new = x.MoveBy(xi.Elements.ToDictionary(kvp => kvp.Key, kvp => t * kvp.Value / xi_norm));
            x_new = x_new.Constrain(area);

            return x_new;
        }
        
        private static RealVector GenerateNewPoint(RealVector x, Area area, double t)
        {
            var xi = GenerateXi(area);
            return MoveByXi(x, xi, area, t);
        }

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
                    var x_new = GenerateNewPoint(x, area, t);
                    
                    var f_x_new = x_new.GetPerformance(f);
                    if (f_x_new < f_x)
                    {
                        var z = x.MoveBy(((x_new - x) * alpha).Elements);
                        var f_z = z.GetPerformance(f);

                        x = x_new;
                        f_x = f_x_new;

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
        
        private class BestSampleNode : GeneralNode<RealVector, double, RealVector>
        {

            public BestSampleNode(int nodeId)
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
                
                var R = state.GetParameter<int>(numberOfGeneratedPointsName);
                var N2 = state.GetParameter<int>(numberOfIterationsOfBestSampleName);
                
                var t_min = state.GetParameter<double>(minStepValueName);
                
                var t = 1.0;

                for (int k = 0; k < N2; ++k)
                {
                    RealVector y_m = null;
                    var f_y_m = double.PositiveInfinity;
                    
                    RealVector[] xi_vectors = new RealVector[R];
                    double[] f_y_values = new double[R];

                    for (int i = 0; i < R; ++i)
                    {
                        xi_vectors[i] = GenerateXi(area);

                        var y = MoveByXi(x, xi_vectors[i], area, t);
                        f_y_values[i] = y.GetPerformance(f);

                        if (f_y_values[i] < f_y_m)
                        {
                            y_m = y;
                            f_y_m = f_y_values[i];
                        }
                    }

                    if (f_y_m < f_x)
                    {
                        x = y_m;
                        f_x = f_y_m;

                        t *= alpha;
                        continue;
                    }
                    else
                    {
                        RealVector stat_grad = area.ToDictionary(kvp => kvp.Key, kvp => 0.0);
                        for (int i = 0; i < R; ++i)
                        {
                            stat_grad = (stat_grad + xi_vectors[i] * (f_y_values[i] - f_x)).Elements;
                        }

                        stat_grad = (stat_grad * (-1.0 / t)).Elements;

                        RealVector y_m_1 = MoveByXi(x, stat_grad, area, t);
                        var f_y_m_1 = y_m_1.GetPerformance(f);

                        if (f_y_m_1 < f_y_m)
                        {
                            x = y_m_1;
                            f_x = f_y_m_1;
                            
                            t *= alpha;
                            continue;
                        }
                    }
                    
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