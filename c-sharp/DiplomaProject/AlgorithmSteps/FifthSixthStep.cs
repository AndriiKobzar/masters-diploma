using System;
using System.Linq;
using System.Threading.Tasks;
using Accord.Math.Integration;

namespace DiplomaProject.AlgorithmSteps
{
    public static class FifthSixthStep
    {
        public static double CoolDelta(double t, int N, double h, double[] sArray, double alpha, double[] y,
            double[] sbmIncrements,
            Func<double, double> sigma,
            Func<double, double> sigmaDerivative,
            Func<double, double> secondSigmaDerivative)
        {
            double eta = 1 / sArray.Sum(x => x * x);
            Func<double, double, double, double> innerIntegrand = Integrand(t, N, h, alpha);
            Func<double, double> steppedY = StepFunction.GetStepFunction(y, t / N);
            var first = 2 * eta * ThirdStep.C(h) * Enumerable.Range(1, N-1).AsParallel().Sum(k => Math.Pow(k * t / N, 0.5 - h) *
                            NonAdaptiveGaussKronrod.Integrate(s => sigma(steppedY(s)) * sigmaDerivative(steppedY(s)) *
                                                                   NonAdaptiveGaussKronrod.Integrate(v => innerIntegrand(s, v, k), k*t/N, s)
                                , k*t/N, t) * sbmIncrements[k]);
            
            var stepSArray = StepFunction.GetStepFunction(sArray, t / N);
            var sumSnDsn = NonAdaptiveGaussKronrod.Integrate(tau => stepSArray(tau) * dvbEta(tau), 0, t);
            return first - sumSnDsn;

            double dvbEta(double v)
            {
                Func<double, double, double> dvbYt = ThirdStep.DvbYt(h, alpha);
                return -4 * Math.Pow(eta, 2) * NonAdaptiveGaussKronrod.Integrate(
                           q => stepSArray(q) * NonAdaptiveGaussKronrod.Integrate(
                                    tau => (Math.Pow(sigmaDerivative(steppedY(tau)), 2) +
                                            sigma(steppedY(tau)) * secondSigmaDerivative(steppedY(tau))) *
                                           dvbYt(v, tau) * dvbYt(q, tau), 0, t), 0, t);
            }
        }

        private static Func<double, double, double, double> Integrand(double t, int n, double h, double alpha)
        {
            return (s, v, k) =>
            {
                double pow = Math.Exp(-alpha * (s - v))
                             * Math.Pow(v, h - 0.5)
                             * Math.Pow(v - k * t / n, h - 1.5);
                return pow;
            };
        }
    }
}