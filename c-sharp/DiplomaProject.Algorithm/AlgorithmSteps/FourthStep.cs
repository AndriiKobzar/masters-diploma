using System;
using System.Linq;
using Accord.Math.Integration;

namespace DiplomaProject.Algorithm.AlgorithmSteps
{
  public static class FourthStep
  {
    public static double[] S(double t, int n, FunctionWithDerivatives sigma, double[] observations, double alpha, double h)
    {
      Func<double, double, double> integrand = DvbSigmaSquareY(t, n, sigma, observations, alpha, h);
      double GetElementOfResult(int k) =>
        2 * IntegrationWrapper.Integrate(x => integrand(k, x), 0,t, "S array integral");

      return Enumerable.Range(0, n)
        .Select(GetElementOfResult)
        .ToArray();
    }

    private static Func<double, double, double> DvbSigmaSquareY(double t, int n, FunctionWithDerivatives sigma,
      double[] observations, double alpha, double h)
    {
      Func<double, double> observationsStep = StepFunction.GetStepFunction(observations, t / n);
      Func<double, double, double> dvbyt = ThirdStep.DvbYt(h, alpha);
      return (k, s) =>
      {
        return sigma.Function(observationsStep(s)) * sigma.FirstDerivative(observationsStep(s)) * dvbyt(k*t/n,s);
      };
    }
  }
}