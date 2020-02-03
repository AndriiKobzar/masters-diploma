using System;
using System.Linq;
using System.Numerics;
using Accord.Math;

namespace DiplomaProject.Algorithm.AlgorithmSteps
{
  public class FgnFeneration
  {
    public static double[] GetFractionalBrownianMotionIncrements(Random random, double t, int n, double h)
    {
      n = (int)Math.Pow(2, Math.Ceiling(Math.Log(n, 2))) + 1;
      var lambda = Lambda(h, n);
      var m = lambda.Length;
      var inverseTransformOnRandom =
        Enumerable.Range(0, m).Select(_ => new Complex(StandardWienerProcess.GaussDistribution(random, 1), 0)).ToArray();
      FourierTransform.FFT(inverseTransformOnRandom, FourierTransform.Direction.Backward);
      var a = lambda.Zip(inverseTransformOnRandom, (x, y) => x * y).ToArray();
      FourierTransform.FFT(a, FourierTransform.Direction.Forward);
      var result = a.Take(m / 2).Select(x => x.Real * Math.Pow(t / n, h)).ToList();
      var average = result.Average();
      return result.Select(x => x - average).ToArray();
    }

    public static double[] Lambda(double h, int n)
    {
      var m = 2 * n - 2;
      var c = new double[m];
      var g = 2 * h;
      double fbc(int x) => (Math.Pow(x + 1, g) + Math.Pow(Math.Abs(x - 1), g) - 2 * Math.Pow(x, g)) / 2;
      for (int i = 0; i < n; i++)
      {
        c[i] = fbc(i);
      }
      var reverseFirstPart = c.Skip(1).Take(n - 2).Reverse().ToArray();
      for (int i = n, j = 0; i < m; i++, j++)
      {
        c[i] = reverseFirstPart[j];
      }
      var result = c.Select(x => new Complex(x, 0)).ToArray();

      FourierTransform.FFT(result, FourierTransform.Direction.Forward);
      return result.Select(x => Math.Sqrt(x.Real)).ToArray();
    }
  }
}
