using System;
using System.Linq;

namespace DiplomaProject.AlgorithmSteps
{
    public static class StandardWienerProcess
    {
        public static double[] GetStandardWienerProcess(Random random, int n, double d)
        {
            return Enumerable.Range(0, n).Select(_ => GaussDistribution(random, d)).ToArray();
        }

        public static double GaussDistribution(Random random, double d)
        {
            int m = 0;
            double x, y, s;
            do
            {
                // випадкові величини рівномірно розподілені на [-1;1]
                x = 2 * random.NextDouble() - 1;
                y = 2 * random.NextDouble() - 1;
                s = x * x + y * y;
            } while (s < 0 || s > 1);

            //  повертає нормально розподілену величину
            return m + Math.Sqrt(d) * x * Math.Sqrt(-2 * Math.Log(s, Math.E) / s);
        }
        
    }
}