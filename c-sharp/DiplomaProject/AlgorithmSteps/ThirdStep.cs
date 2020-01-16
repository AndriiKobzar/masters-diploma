using System;
using System.Linq;
using System.Threading.Tasks;
using Accord.Math.Integration;

namespace DiplomaProject.AlgorithmSteps
{
    public static class ThirdStep
    {
        public static Func<double, double, double> DvbYt(double h, double alpha)
        {
            Func<double, double, double> kernel = Kernel(h);
            return (u, t) => kernel(t, u) - alpha * Math.Exp(-alpha * t) *
                             NonAdaptiveGaussKronrod.Integrate(s => Math.Exp(alpha * s) * kernel(s, u), u, t);
        }

        private static Func<double, double, double> Kernel(double h)
        {
            return (s, t) => s < t
                ? C(h) * Math.Pow(s, 0.5 - h) *
                  NonAdaptiveGaussKronrod.Integrate(
                      u => Math.Pow(u, h - 0.5) * Math.Pow(u - s, h - 1.5),
                      s, t)
                : 0;
        }


        public static double C(double h)
        {
            double dividend = 2 * h * Gamma(1.5 - h);
            double divisor = Gamma(h + 0.5) * Gamma(2 - 2 / h);
            return (h - 0.5) * Math.Sqrt(dividend / divisor);
        }

        private static double Gamma(double x) => Accord.Math.Gamma.Function(x);
    }
}