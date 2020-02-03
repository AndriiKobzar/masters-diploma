using System;
using System.Linq;

namespace DiplomaProject.Algorithm.AlgorithmSteps
{
    public class Density
    {
        public static double GetDensity(Random random, double h, double t, int n, double alpha,
            FunctionWithDerivatives sigma, double u)
        {
            double[] fbmIncrements = EulerScheme.GetFractionalBrownianMotionIncrements(random, h, n, t);
            double[] observations = EulerScheme.GetObservations(alpha, n, t, fbmIncrements).ToArray();
            double integral = SecondStep.IntegralOfSquareSigma(t, sigma.Function, observations);

            if (integral <= u)
            {
                return 0;
            }

            double[] s = FourthStep.S(t, n, sigma, observations, alpha, h);
            double[] sbmIncrements = StandardWienerProcess.GetStandardWienerProcess(random, n, t / n);
            double delta = FifthSixthStep.CoolDelta(t, n, h, s, alpha, observations, sbmIncrements,
                sigma);
            return delta;
        }
    }
}