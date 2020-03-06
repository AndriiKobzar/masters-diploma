using System;
using System.Linq;
using Accord.Math.Integration;

namespace DiplomaProject.Algorithm.AlgorithmSteps
{
    public static class FifthSixthStep
    {
        public static double CoolDelta(double t, int N, double h, double[] sArray, double alpha, double[] y,
            double[] sbmIncrements, FunctionWithDerivatives sigma)
        {
            double eta = 1 / sArray.Sum(x => x * x);
            Func<double, double, double, double> innerIntegrand = Integrand(t, N, h, alpha);
            Func<double, double> steppedY = StepFunction.GetStepFunction(y, t / N);
            var first = 2 * eta * ThirdStep.C(h) * Enumerable.Range(1, N-1).AsParallel().Sum(k => Math.Pow(k * t / N, 0.5 - h) *
                            IntegrationWrapper.Integrate(s => sigma.Function(steppedY(s)) * sigma.FirstDerivative(steppedY(s)) *
                                                                   IntegrationWrapper.Integrate(v => innerIntegrand(s, v, k), 0, s, "delta first dv")
                                , 0, t, "delta first ds") * sbmIncrements[k]);
            
            var stepSArray = StepFunction.GetStepFunction(sArray, t / N);
            var sumSnDsn = IntegrationWrapper.Integrate(tau => stepSArray(tau) * dvbEta(tau), 0, t, "sn*dsn");
            return first - sumSnDsn;

            double dvbEta(double v)
            {
                Func<double, double, double> dvbYt = ThirdStep.DvbYt(h, alpha);
                return -4 * Math.Pow(eta, 2) * IntegrationWrapper.Integrate(
                           q => stepSArray(q) * IntegrationWrapper.Integrate(
                                    tau => (Math.Pow(sigma.FirstDerivative(steppedY(tau)), 2) +
                                            sigma.Function(steppedY(tau)) * sigma.SecondDerivative(steppedY(tau))) *
                                           dvbYt(v, tau) * dvbYt(q, tau), 0, t, "dvbeta dtau"), 0, t, "dvbeta dq");
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