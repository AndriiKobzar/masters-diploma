using System;
using Accord.Math.Integration;

namespace DiplomaProject
{
    public class PriceForecast
    {
        static double Forecast(double[] historicData, double[] density, double t, double b, double n, Func<double, double> payoffFunction)
        {
            Func<double, double> stepDensity = StepFunction.GetStepFunction(density, t/n);
            return 2 * Math.Pow(Math.PI, -0.5)*t*NonAdaptiveGaussKronrod.Integrate(x => G(x) * 
                       NonAdaptiveGaussKronrod.Integrate(u => (x+u/2-historicData[0] - b*t)/Math.Pow(u,3) *
                                                              Math.Pow(Math.Exp(-(x+u/2-historicData[0] - b*t)),2)/(2*Math.Pow(u,2))*stepDensity(u), 0,t), 0, t);


            double G(double x)
            {
                return NonAdaptiveGaussKronrod.Integrate(z => payoffFunction(Math.Exp(z)), 0, x);
            }
        }
        
    }
}