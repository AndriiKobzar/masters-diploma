using System;
using System.Linq;

namespace DiplomaProject.AlgorithmSteps
{
    public static class SecondStep
    {
        /// <summary>
        /// Calculates Integral of square sigma for observations array (second step)
        /// </summary>
        /// <param name="t">Length</param>
        /// <param name="sigma">Sigma square function</param>
        /// <param name="observations">Array of observations</param>
        /// <returns></returns>
        public static double IntegralOfSquareSigma(double t, Func<double, double> sigma, double[] observations)
        {
            double delta = t / observations.Length;
            return observations.Select(y => Math.Pow(sigma(y), 2) * delta).Sum();
        } 
    }
}