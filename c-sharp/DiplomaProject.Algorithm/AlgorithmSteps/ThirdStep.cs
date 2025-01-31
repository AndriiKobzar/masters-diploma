using System;
using Accord.Math.Integration;

namespace DiplomaProject.Algorithm.AlgorithmSteps
{
    public static class ThirdStep
    {
        public static Func<double, double, double> DvbYt(double h, double alpha)
        {
            // Func<double, double, double> kernel = Kernel(h);
            // return (u, t) => kernel(t, u) - alpha * Math.Exp(-alpha * t) *
            //                  IntegrationWrapper.Integrate(s => Math.Exp(alpha * s) * kernel(s, u), u, t);
            return (u, t) =>
            {
                return u < t
                    ? 
                      IntegrationWrapper.Integrate(
                          s => C(h) * Math.Exp(-alpha * t) * Math.Pow(u, 0.5 - h) * Math.Exp(alpha * s) * Math.Pow(s, h - 0.5) * Math.Pow(s - u, h - 1.5), u, t, "dvbyt")
                    : 0;
            };
        }

        // private static Func<double, double, double> Kernel(double h)
        // {
        //     return (s, t) => s < t
        //         ? C(h) * Math.Pow(s, 0.5 - h) *
        //           IntegrationWrapper.Integrate(
        //               u => Math.Pow(u, h - 0.5) * Math.Pow(u - s, h - 1.5),
        //               s, t)
        //         : 0;
        // }


        public static double C(double h)
        {
            double dividend = 2 * h * Gamma(1.5 - h);
            double divisor = Gamma(h + 0.5) * Gamma(2 - 2 / h);
            return (h - 0.5) * Math.Sqrt(dividend / divisor);
        }

        private static double Gamma(double x) => Accord.Math.Gamma.Function(x);
    }
}