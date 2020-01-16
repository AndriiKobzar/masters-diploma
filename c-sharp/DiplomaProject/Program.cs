using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using DiplomaProject.AlgorithmSteps;

namespace DiplomaProject
{
    class Program
    {
        static void Main()
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


            var pointsCount = 50;
            var rangeStart = 0;
            var rangeEnd = 2;
            ConcurrentDictionary<double, double> resultSet = new ConcurrentDictionary<double, double>();
            var fileName = string.Format("diploma-log-{0}.txt", DateTime.Now);
            //using (var streamWriter = new StreamWriter(new FileStream(fileName, FileMode.OpenOrCreate)))
            {
                //Parallel.For(1, pointsCount+1, (i) =>
                //{
                for (int i = 1; i <= pointsCount; i++)
                {
                    var random = new Random(DateTime.Now.Millisecond);
                    double u = rangeStart + 1.0 * i * rangeEnd / pointsCount;
                    double averageDelta = Enumerable.Range(1, 100)
                      .AsParallel()
                      .WithDegreeOfParallelism(100)
                      .Select(calculationIndex =>
                      {
                          double[] observations;
                          double[] fbmIncrements;
                          double integral;

                          //do
                          //{
                              Console.WriteLine("Calculating Y array {0} {1}", u, calculationIndex);
                              fbmIncrements = EulerScheme.GetFractionalBrownianMotionIncrements(random, h, n, t);
                              observations = EulerScheme.GetObservations(alpha, n, t, fbmIncrements).ToArray();
                              integral = SecondStep.IntegralOfSquareSigma(t, Sigma, observations);
                              Console.WriteLine("integral of sigma {0} {1} : {2}", u, calculationIndex, integral);
                          //} while (integral > u);

                          if (integral <= u)
                          {
                              Console.WriteLine("delta {0} {1} = {2}", u, calculationIndex, 0);
                              return 0;
                          }

                          //Console.WriteLine("fbm increments");
                          //Display.Array(fbmIncrements);

                          //Console.WriteLine("Y array");
                          //Display.Array(observations);

                          Console.WriteLine("Calculating array S {0} {1}", u, calculationIndex);
                          double[] s = FourthStep.S(t, n, Sigma, SigmaDerivative, observations, alpha, h);
                          Console.WriteLine("Calculated array S {0} {1}", u, calculationIndex);
                          //Display.Array(s);
                          Console.WriteLine("Calculating delta {0} {1}", u, calculationIndex);
                          double[] sbmIncrements = StandardWienerProcess.GetStandardWienerProcess(random, n, t / n);
                          double delta = FifthSixthStep.CoolDelta(t, n, h, s, alpha, observations, sbmIncrements,
                    Sigma,
                    SigmaDerivative,
                    SigmaSecondDerivative);
                          Console.WriteLine("delta {0} {1} = {2}", u, calculationIndex, delta);
                         // streamWriter.WriteLine("delta {0} {1} = {2}", u, calculationIndex, delta);
                          return delta;
                      }).Average();
                    resultSet.TryAdd(u, averageDelta);
                    Console.WriteLine(DateTime.Now);
                }
                Console.WriteLine(DateTime.Now);
            }
            //});
            //Console.ReadLine();
            //await new WolframClient().ExecuteRequest("D[sin[x]^2, {x, 2}]");
            foreach (KeyValuePair<double, double> pair in resultSet)
            {
                Console.WriteLine($"{pair.Key} {pair.Value}");
            }
        }
    }
}