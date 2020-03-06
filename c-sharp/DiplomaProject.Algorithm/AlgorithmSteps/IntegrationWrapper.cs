using System;
using Accord;
using Accord.Math.Integration;

namespace DiplomaProject.Algorithm.AlgorithmSteps
{
    public class IntegrationWrapper
    {
        public static double Integrate(Func<double, double> func, double a, double b, string name)
        {
            
            var integrator = new InfiniteAdaptiveGaussKronrod(12000000)
            {
                Function = func,
                Range = new DoubleRange(a, b)
            };
            integrator.Compute();
            if (integrator.Status != InfiniteAdaptiveGaussKronrodStatus.Success)
            {
                //Console.WriteLine($"{name} {integrator.Status} result {integrator.Area} error {integrator.Error}. Boundaries: {a} {b}");
            }
            return integrator.Area;
        }
    }
}