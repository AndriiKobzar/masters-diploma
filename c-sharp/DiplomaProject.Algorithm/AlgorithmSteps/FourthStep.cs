using System;
using System.Linq;
using Accord.Math.Integration;

namespace DiplomaProject.Algorithm.AlgorithmSteps
{
  public static class FourthStep
  {
    public static double[] S(double t, int N, FunctionWithDerivatives sigma, double[] observations, double alpha, double h)
    {
      Func<double, double, double> integrand = DvbSigmaSquareY(t, N, sigma, observations, alpha, h);
      double GetElementOfResult(int k) =>
        2 * NonAdaptiveGaussKronrod.Integrate(x => integrand(k, x), 0,t);

      return Enumerable.Range(0, N)
        .AsParallel()
        .WithDegreeOfParallelism(4)
        .Select(GetElementOfResult)
        .ToArray();
    }

    public static Func<double, double, double> DvbSigmaSquareY(double t, int N, FunctionWithDerivatives sigma,
      double[] observations, double alpha, double h)
    {
      Func<double, double> observationsStep = StepFunction.GetStepFunction(observations, t / N);
      return (k, s) =>
      {
        Func<double, double, double> dvbyt = ThirdStep.DvbYt(h, alpha);
        return sigma.Function(observationsStep(s)) * sigma.FirstDerivative(observationsStep(s)) * dvbyt(k*t/N,s);
      };
    }
  }
}