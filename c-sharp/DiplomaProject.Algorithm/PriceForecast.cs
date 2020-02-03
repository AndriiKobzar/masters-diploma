using System;
using System.Linq;
using Accord.Math.Integration;

namespace DiplomaProject.Algorithm
{
    public class PriceForecast
    {
        public static double Forecast(double firstHistoricValue, Point[] density, double b, Func<double, double> payoffFunction, double t)
        {
            var maxDensity = density.Select(p => p.X).Max();
            Func<double, double> stepDensity = StepFunction.GetStepFunction(density);
            return 2 * Math.Pow(Math.PI, -0.5)*t*NonAdaptiveGaussKronrod.Integrate(x => G(x) * 
                       NonAdaptiveGaussKronrod.Integrate(u => (x+u/2-firstHistoricValue - b*t)/Math.Pow(u,3) *
                                                              Math.Pow(Math.Exp(-(x+u/2-firstHistoricValue - b*t)),2)/(2*Math.Pow(u,2))*stepDensity(u), 0,maxDensity), 0, t);


            double G(double x)
            {
                return NonAdaptiveGaussKronrod.Integrate(z => payoffFunction(Math.Exp(z)), 0, x);
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