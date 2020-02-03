using System;

namespace DiplomaProject.Algorithm
{
    public class FunctionWithDerivatives
    {
        public Func<double, double> Function { get; set; }
        public Func<double, double> FirstDerivative { get; set; }
        public Func<double, double> SecondDerivative { get; set; }
    }
}