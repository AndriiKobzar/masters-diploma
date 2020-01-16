using System;
using System.Numerics;

namespace DiplomaProject.AlgorithmSteps
{
  public class fG
  {
    public double H { get; set; }
    public int Q { get; set; }
    public int N { get; set; }
    public int n { get; set; }
    public Complex[] C;
    public Complex[] fft(Complex[] input)
    {
      Complex[] res = new Complex[input.Length];
      Complex e = 0;
      double koeff = 0;
      Complex sum = 0;
      for (int k = 0; k < input.Length; k++)
      {
        for (int u = 0; u < input.Length; u++)
        {
          koeff = -2 * Math.PI * u * k / input.Length;

          e = new Complex(Math.Cos(koeff), Math.Sin(koeff));
          sum += (input[u] * e);
        }
        res[k] = sum;
        sum = 0;
      }
      return res;
    }
    public Complex[] ifft(Complex[] input)
    {
      Complex[] res = new Complex[input.Length];
      Complex e = 0;
      double koeff = 0;
      Complex sum = 0;
      for (int k = 0; k < input.Length; k++)
      {
        for (int u = 0; u < input.Length; u++)
        {
          koeff = 2 * Math.PI * u * k / input.Length;

          e = new Complex(Math.Cos(koeff), Math.Sin(koeff));
          sum += (input[u] * e);
        }
        res[k] = sum / input.Length;
        sum = 0;
      }
      return res;
    }
    public Complex[] fbc(Complex[] res, double G, int N)
    {
      for (int j = 0; j < N; j++)
      {
        res[j] = (Math.Pow((j + 1), G) + Math.Pow(Math.Abs(j - 1), G) - 2 * Math.Pow(j, G)) / 2;
      }
      return res;
    }
    public Complex[] fliplr(Complex[] input)
    {
      Complex temp = 0;
      for (int i = 0; i < Math.Floor(input.Length / 2.0); i++)
      {
        temp = input[i];
        input[i] = input[input.Length - i - 1];
        input[input.Length - i - 1] = temp;
      }
      return input;
    }

    public Complex[] CreateVector(double H, int N)
    {
      int M = 2 * N - 2;
      Complex[] C = new Complex[M];
      return C;

    }
    public Complex[] Swap(Complex[] input, int N, Complex[] littleMass)
    {
      int i = 0;
      Complex[] res = input;
      for (int j = N; j < 2 * N - 2; j++)
      {
        res[j] = littleMass[i];
        i++;
      }
      return res;
    }
    public Complex[] Part(Complex[] input, int N)
    {
      int lenght = N - 2;
      Complex[] res = new Complex[lenght];
      int j = 0;
      for (int i = 1; i < N - 1; i++)
      {
        res[j] = input[i];
        j++;
      }
      return res;
    }
    public int Size(Complex[] input)
    {
      int res = input.Length;
      return res;
    }
    public Complex[,] CreateMatrix(Random rand, int k, int l)
    {
      int sign = 0;
      Complex[,] res = new Complex[k, l];
      for (int i = 0; i < k; i++)
      {
        for (int j = 0; j < l; j++)
        {
          sign = rand.Next(-2, 2);
          if (sign == 0)
          {
            res[i, j] = rand.NextDouble();
          }
          else
          {
            res[i, j] = rand.NextDouble() * sign;
          }
        }
      }
      return res;
    }
    public Complex[] PowMatrix(Complex[] input, double pow)
    {
      for (int i = 0; i < input.Length; i++)
      {
        input[i] = Math.Pow(input[i].Real, pow);
      }
      return input;
    }
    public Complex[,] MultiplyMatrix(Complex[,] firstMatrix, Complex[] secondMatrix)
    {
      for (int i = 0; i < firstMatrix.Length / secondMatrix.Length; i++)
      {
        for (int j = 0; j < secondMatrix.Length; j++)
        {
          firstMatrix[i, j] = firstMatrix[i, j] * secondMatrix[j];
        }
      }
      return firstMatrix;
    }
    public double[,] TakePartOfMatrix(Complex[,] input, int M)
    {
      double[,] res = new double[input.Length / M, M / 2];

      for (int i = 0; i < input.Length / M; i++)
      {
        for (int j = 0; j < M / 2; j++)
        {
          res[i, j] = input[i, j].Real;
        }
      }
      return res;
    }
    public Complex[,] SetRowIntoMatrix(Complex[,] res, Complex[] input, int count)
    {
      for (int i = 0; i < input.Length; i++)
      {
        res[count, i] = input[i];
      }
      return res;
    }
    public Complex[,] IfftAndFftForMatrix(Complex[,] input, int M, int n, bool flag)
    {
      Complex[,] res = new Complex[n, M];
      Complex[] row = new Complex[M];
      for (int i = 0; i < input.Length / M; i++)
      {
        for (int j = 0; j < M; j++)
        {
          row[j] = input[i, j];
        }
        if (flag)
        {
          res = SetRowIntoMatrix(res, fft(row), i);
        }
        else
        {
          res = SetRowIntoMatrix(res, ifft(row), i);
        }
      }
      return res;
    }

    public Complex[] Lambda(double _h, int _q)
    {
      H = _h;
      Q = _q;
      N = (int)Math.Pow(2, Q) + 1;
      C = CreateVector(H, N);
      Complex[] temp;
      C = fbc(C, 2 * H, N);
      temp = Part(C, N);
      temp = fliplr(temp);
      C = Swap(C, N, temp);
      C = fft(C);
      C = PowMatrix(C, 0.5);
      return C;
    }
    public Complex[,] RandomMatrix(Random random, Complex[] l, int _n)
    {
      n = _n;
      int M = Size(l);
      Complex[,] randomMatrix = CreateMatrix(random,n, M);
      return randomMatrix;
    }
    /// <summary>
    /// Returns n paths for fgn
    /// </summary>
    /// <param name="lambda"></param>
    /// <param name="_n"></param>
    /// <param name="randomMatrix"></param>
    /// <returns></returns>
    public double[,] FGN(Complex[] lambda, int _n, Complex[,] randomMatrix)
    {
      int M = Size(lambda);
      randomMatrix = IfftAndFftForMatrix(randomMatrix, M, _n, false);
      randomMatrix = MultiplyMatrix(randomMatrix, lambda);
      randomMatrix = IfftAndFftForMatrix(randomMatrix, M, _n, true);
      double[,] result = TakePartOfMatrix(randomMatrix, M);
      return result;
    }
    public void ShowLambda(Complex[] input)
    {
      Console.Write("Lambda = ");
      for (int i = 0; i < input.Length; i++)
      {
        Console.Write(Math.Round(input[i].Real, 2) + " ");
      }
      Console.WriteLine();
    }
    public void ShowMatrix(Complex[,] input, int _n)
    {
      double inputNumber;
      int count = 0;
      Console.WriteLine("Ramdom Matrix");
      for (int i = 0; i < _n; i++)
      {
        for (int j = 0; j < input.Length / _n; j++)
        {
          inputNumber = Math.Round(input[i, j].Real, 2);
          count = 5 - inputNumber.ToString().Length;
          if (count == 1 || count == 2)
          {
            Console.Write(" " + inputNumber + " ");
          }
          else if (count == 3 || count == 4)
          {
            Console.Write("  " + inputNumber + " ");
          }
          else if (count == 0)
          {
            Console.Write(inputNumber + " ");
          }
          count = 0;
        }
        Console.WriteLine();
      }
      Console.WriteLine();
    }
    public void ShowResult(double[,] input, int _n)
    {
      double inputNumber;
      int count = 0;
      Console.WriteLine("fGnsamples =");
      for (int i = 0; i < _n; i++)
      {
        for (int j = 0; j < input.Length / _n; j++)
        {
          inputNumber = Math.Round(input[i, j], 2);
          count = 5 - inputNumber.ToString().Length;
          if (count == 1 || count == 2)
          {
            Console.Write(" " + inputNumber + " ");
          }
          else if (count == 3 || count == 4)
          {
            Console.Write("  " + inputNumber + " ");
          }
          else if (count == 0)
          {
            Console.Write(inputNumber + " ");
          }
          count = 0;
        }
        Console.WriteLine();
      }
      Console.WriteLine();
    }
    public Complex[] Init(double h, int q)
    {
      Complex[] lambda;
      if (h <= 0 || h >= 1)
      {
        throw new Exception();
      }
      lambda = Lambda(h, q);
      return lambda;

    }
  }
}
