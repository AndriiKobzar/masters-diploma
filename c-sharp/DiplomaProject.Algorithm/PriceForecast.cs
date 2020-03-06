using System;
using System.Linq;
using Accord.Math.Integration;
using DiplomaProject.Algorithm.AlgorithmSteps;

namespace DiplomaProject.Algorithm
{
    public class PriceForecast
    {
        public static double Forecast(double firstHistoricValue, Point[] density, double b,
            Func<double, double> payoffFunction, double t)
        {
            double minDensity = density.Select(p => p.X).Min();
            double maxDensity = density.Select(p => p.X).Max();
            Func<double, double> stepDensity = StepFunction.GetStepFunction(density);
            return Math.Pow(2 * Math.PI, -0.5) * t *
                   IntegrationWrapper.Integrate(x => G(x) * IntegrationWrapper.Integrate(u =>
                                                                  (x + u / 2 - firstHistoricValue - b * t) /
                                                                  Math.Pow(u, 3) *
                                                                  Math.Exp(
                                                                      -Math.Pow(x + u / 2 - firstHistoricValue - b * t,
                                                                          2) / (2 * Math.Pow(u, 2))) * Math.Abs(stepDensity(u)),
                                                              minDensity, maxDensity, "forecast integral density"), -t, t, "forecast integral dx");


            double G(double x)
            {
                return IntegrationWrapper.Integrate(z => payoffFunction(Math.Exp(z)), 0, x, "integral G(x)");
            }
        }
    }

    public struct Point
    {
        public Point(double x, double y)
        {
            X = x;
            Y = y;
        }

        public double X { get; }
        public double Y { get; }
    }
}