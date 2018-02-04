using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Assignment3
{
  class Program
  {
    public static readonly int numIter = 3;
    static int numClass = 3;
    public static readonly double normalDensityConst = Math.Sqrt(1/(2 * 3.14));
    public static double[] input = new double[] {9,10,11,20,21,22,46,49,55,57}; //{-6,-5,-4,0,4,5,6}; // 
    public static double[][] mean = new double[numIter + 1][]; // {-20,6}; // {10,21,51.75};
    static double[][][] normalDensityMatrix = new double[numIter][][];
    static double[][] nomralDensitySumMatrix = new double[numIter][];
    static double[][][] expectedMatrix = new double[numIter][][];
    static void Main(string[] args)
    {
      initialize();
      Console.WriteLine(logLikelihood(0));
      Console.WriteLine(computeBIC(-116.1755));
      Console.WriteLine("Press any key to exit....");
      Console.ReadLine();
    }

    public static void initialize() {
      mean[0] = new double[] {10,21,51.75};
      numClass = mean[0].Length;

      for (int iter = 0; iter < numIter; iter++) {
       // E-Step
        normalDensityMatrix[iter] = new double[input.Length][];
        nomralDensitySumMatrix[iter] = new double[input.Length];
        for (int row = 0; row < input.Length; row++) {
          normalDensityMatrix[iter][row] = new double[numClass];
          nomralDensitySumMatrix[iter][row] = 0;
          for (int col = 0; col < numClass; col++) {
            normalDensityMatrix[iter][row][col] =  normalDensity(input[row], mean[iter][col]);
            nomralDensitySumMatrix[iter][row] += normalDensityMatrix[iter][row][col];
          }
        }
        expectedMatrix[iter] = new double[input.Length][];
        for (int row = 0; row < input.Length; row++) {
          expectedMatrix[iter][row] = new double[numClass];  
          for (int col = 0; col < numClass; col++) {
            expectedMatrix[iter][row][col] =  normalDensityMatrix[iter][row][col]/nomralDensitySumMatrix[iter][row];
          }
        }

        // M-Step
        int idx = iter + 1;
        mean[idx] = new double[numClass];
        for (int col = 0; col < numClass; col++) {
          double denom = 0;
          for (int row = 0; row < input.Length; row++) {
            mean[idx][col] += expectedMatrix[iter][row][col] * input[row];
            denom += expectedMatrix[iter][row][col];
          }
          mean[idx][col] /= denom;
        }
      }
    }

    public static double logLikelihood(int iterNum, double stddev = 1.0) {
      double sum = 0.0;
      for (int row = 0; row < input.Length; row++) {
        sum += Math.Log(1.0/numClass) - 0.5*Math.Log(2 * 3.14 * Math.Pow(stddev,2));
        for(int col =0; col < numClass; col++) {
          sum -= (expectedMatrix[iterNum][row][col]* (input[row] - mean[iterNum + 1][col])/(2 * Math.Pow(stddev,2)));
        }
      }
      return sum;
    }

    public static double computeBIC(double logLik) {
      return 2 * logLik - numClass * Math.Log(input.Length);
    }

    public static double normalDensity(double x, double mean, double stddev = 1)
    {
      double diffSqr = Math.Pow((double)(x - mean), 2.0);
      double stddevSqr = 2 * Math.Pow((double)stddev, 2.0);
      double power = Math.Exp(-1.0 * diffSqr / stddevSqr);
      return (normalDensityConst * power);
    }

  }
}
