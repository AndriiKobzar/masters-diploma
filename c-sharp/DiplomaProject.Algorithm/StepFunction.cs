using System;
using System.Linq;

namespace DiplomaProject.Algorithm
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
        public static Func<double, double> GetStepFunction(Point[] points)
        {
            return x =>
            {
                (double _, Point point) = points
                    .Select(p => (p.X - x, p))
                    .Where(tuple => tuple.Item1 < 0)
                    .OrderBy(tuple => tuple.Item1)
                    .LastOrDefault();
                return point.Y;
            };
        }
    }
}