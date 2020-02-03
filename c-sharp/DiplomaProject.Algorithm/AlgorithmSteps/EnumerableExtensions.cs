using System;
using System.Collections.Generic;
using System.Linq;

namespace DiplomaProject.Algorithm.AlgorithmSteps
{
  public static class EnumerableHelper
  {
    internal static IEnumerable<int> Interval(int start, int end)
    {
      if(start > end)
      {
        throw new ArgumentException("Start should be greater then or equal end");
      }
      return Enumerable.Range(start, end - start + 1);
    }
  }
}
