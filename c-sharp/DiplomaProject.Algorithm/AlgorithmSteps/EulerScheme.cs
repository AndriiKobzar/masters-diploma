using System;
using System.Collections.Generic;
using System.Linq;

namespace DiplomaProject.AlgorithmSteps
{
    public static class EulerScheme
    {
        /// <summary>
        /// Implements Euler scheme for generating array of observations
        /// </summary>
        /// <param name="alpha">Alpha parameter</param>
        /// <param name="n"></param>
        /// <param name="t"></param>
        /// <param name="fbmIncrements"></param>
        /// <returns></returns>
        public static IEnumerable<double> GetObservations(double alpha, int n, double t,
            double[] fbmIncrements)
        {
            double incrementOfN = t / n;
            double currentValue = 0;
            yield return 0;
            for (int i = 0; i < n; i++)
            {
                currentValue = currentValue - alpha * currentValue * incrementOfN + fbmIncrements[i];
                
                yield return currentValue;
            }
        }

        public static double[] GetFractionalBrownianMotionIncrements(Random random, double h, int n, double t)
        {
            //return FgnFeneration.GetFractionalBrownianMotionIncrements(random, t, n, h);
            var q = (int) Math.Ceiling(Math.Log(n - 1, 2));
            return GetFractionalGaussianNoise(random, h, q).Select(s => s * Math.Pow(t / n, h)).ToArray();
        }

        private static IEnumerable<double> GetFractionalGaussianNoise(Random rand, double h, int q)
        {
            var fg = new fG();
            var lambda = fg.Init(h, q);

            var fractionalGaussianNoise = fg.FGN(lambda, 1, fg.RandomMatrix(rand, lambda, 1));

            for (int i = 0; i < fractionalGaussianNoise.GetLength(1); i++)
            {
                yield return fractionalGaussianNoise[0, i];
            }
        }
    }
}