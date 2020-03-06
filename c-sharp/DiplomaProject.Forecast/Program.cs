using System;
using System.IO;
using System.Linq;
using DiplomaProject.Algorithm;

namespace DiplomaProject.Forecast
{
    class Program
    {
        static void Main(string[] _)
        {
            Point[] density = GetDensityData("density.txt");
            double x0 = 0;
            Console.WriteLine(0);
            const double b = 0.4;
            const double t = 1;
            static double PayoffFunction(double s) => Math.Max(s - 1, 0) + (s > 1 ? 1 : 0);
            for (int i = 0; i < 100; i++)
            {
                x0 = PriceForecast.Forecast(x0, density, b, PayoffFunction, t);
                Console.WriteLine(x0);
            }
        }

        static Point[] GetDensityData(string fileName)
        {
            return File.ReadLines(fileName)
                .Select(line => line.Split(' ').Select(double.Parse).ToArray())
                .Select(numbers => new Point(numbers[0], numbers[1]))
                .ToArray();
        }
    }
}