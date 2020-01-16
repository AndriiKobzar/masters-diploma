using System;

namespace DiplomaProject
{
    public static class StepFunction
    {
        public static Func<double, double> GetStepFunction(double[] points, double interval)
        {
            return x =>
            {
                int k = (int) (x / interval);
                return points[k];
            };
        }
    }
}