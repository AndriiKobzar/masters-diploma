using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using DiplomaProject.Algorithm;
using DiplomaProject.Algorithm.AlgorithmSteps;

namespace DiplomaProject
{
    public class Program
    {
        public static void Main()
        {
            Console.WriteLine(DateTime.Now);
            const int n = 250;
            const double t = 1.0;
            const double h = 0.6;
            const double alpha = 0.6;

            double Sigma(double y) => Math.Pow(Math.Sin(y), 2) + 0.05;
            double SigmaDerivative(double y) => Math.Sin(2 * y);
            double SigmaSecondDerivative(double y) => 2 * Math.Cos(2 * y);

            //double Sigma(double y) => Math.Sqrt(y * y + 1);
            //double SigmaDerivative(double y) => y * Math.Pow(y * y + 1, -0.5);
            //double SigmaSecondDerivative(double y) => Math.Pow(y * y + 1, -0.5) - y * y * Math.Pow(y * y + 1, -1.5);


            double PayoffFunction(double s) => Math.Max(s - 1, 0) + (s > 1 ? 1 : 0);

            var pointsCount = 50;
            var rangeStart = 0;
            var rangeEnd = 2;
            ConcurrentDictionary<double, double> resultSet = new ConcurrentDictionary<double, double>();
            var sigma = new FunctionWithDerivatives
            {
                Function = Sigma,
                FirstDerivative = SigmaDerivative,
                SecondDerivative = SigmaSecondDerivative
            };
            for (int i = 1; i <= pointsCount; i++)
            {
                var random = new Random(DateTime.Now.Millisecond);
                double u = rangeStart + 1.0 * i * rangeEnd / pointsCount;
                double averageDelta = Enumerable.Range(1, 100)
                    .AsParallel()
                    .WithDegreeOfParallelism(100)
                    .Select(calculationIndex => Density.GetDensity(random, h, t, n, alpha, sigma, u))
                    .Average();
                resultSet.TryAdd(u, averageDelta);
                Console.WriteLine(DateTime.Now);
            }

            Console.WriteLine(DateTime.Now);
            
            foreach (KeyValuePair<double, double> pair in resultSet)
            {
                Console.WriteLine($"{pair.Key} {pair.Value}");
            }
        }
    }
}