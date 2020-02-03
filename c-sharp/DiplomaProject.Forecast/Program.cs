using System;
using System.IO;
using System.Linq;
using DiplomaProject.Algorithm;

namespace DiplomaProject.Forecast
{
    class Program
    {
        static void Main(string[] args)
        {
            Point[] density = GetDensityData("density.txt");
            const int x0 = 0;
            const double b = 0.2;
            static double PayoffFunction(double s) => Math.Max(s - 1, 0) + (s > 1 ? 1 : 0);
            double forecast = PriceForecast.Forecast(x0,density,b,PayoffFunction,1);
            Console.WriteLine(forecast);
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