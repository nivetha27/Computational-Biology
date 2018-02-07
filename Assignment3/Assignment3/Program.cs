using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Assignment3
{
  class Program
  {
    public static readonly int numIter = 11;
    static int numClass = 5;
    public static readonly double normalDensityConst = Math.Sqrt(1/(2 * Math.PI));
    public static double[] input = new double[]  {9,10,11,20,21,22,46,49,55,57};  // {-6,-5,-4,0,4,5,6}; 
    public static double[][] meanArr = new double[numIter + 1][];
    public static double[] stddevArr = new double[numClass];
    static double[][][] normalDensityMatrix = new double[numIter][][];
    static double[][] nomralDensitySumMatrix = new double[numIter][];
    static double[][][] expectedMatrix = new double[numIter][][];
    static void Main(string[] args)
    {      
      initializeCluster();
      initialize();
      // Console.WriteLine("Press any key to exit....");
      // Console.ReadLine();
      print(meanArr);
    }

    public static void initializeCluster() {
      //int numSamples = input.Length/numClass;
      //meanArr[0] = new double[2] { -20.0, 6.0}; // new double[numClass];
      //for (int i = 0; i < numClass; i++) {
      //  int start = i * numSamples;
      //  int end = (start + numSamples);
      //  if (i == numClass - 1) {
      //    end = input.Length;  
      //  }
      //  // meanArr[0][i] =  computeMean(start, end);
      //  stddevArr[i] = 1.0; //computeStddev(start, end, i);
      //}

      //var rand = new Random();
      //  meanArr[0] = new double[numClass];
      //var min = input.Min();
      //var max = input.Max();

      //for (int i = 0; i < numClass; i++)
      //{
      //    meanArr[0][i] = rand.Next((int)min, (int)max);
      //    stddevArr[i] = 1.0;
      //}
      meanArr[0] = new double[5] { 35,12,46,22,45  }; 
      stddevArr = new double[5] {1,1,1,1,1};
    }

    public static void initialize() {
      for (int iter = 0; iter < numIter; iter++) {
       // E-Step
        normalDensityMatrix[iter] = new double[input.Length][];
        nomralDensitySumMatrix[iter] = new double[input.Length];
        for (int row = 0; row < input.Length; row++) {
          normalDensityMatrix[iter][row] = new double[numClass];
          nomralDensitySumMatrix[iter][row] = 0;
          for (int col = 0; col < numClass; col++) {
            normalDensityMatrix[iter][row][col] =  normalDensity(input[row], meanArr[iter][col],stddevArr[col]);
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
        meanArr[idx] = new double[numClass];
        for (int col = 0; col < numClass; col++) {
          double denom = 0;
          for (int row = 0; row < input.Length; row++) {
            meanArr[idx][col] += expectedMatrix[iter][row][col] * input[row];
            denom += expectedMatrix[iter][row][col];
          }
          meanArr[idx][col] /= denom;
        }

        double logLikelihood = (computeLogLikelihood(iter));
        Console.Write(logLikelihood);
        Console.WriteLine(computeBIC(logLikelihood));
      }
    }

    public static double computeLogLikelihood(int iterNum, double stddev = 1.0) {
      double logLikelihood = 0;
      for (int i = 0; i < normalDensityMatrix[0].Length; i++)
      {
        double sum = 0;
        for (int j = 0; j < normalDensityMatrix[0][0].Length; j++)
        {
          sum += (1.0 / numClass * normalDensityMatrix[iterNum][i][j]);
        }
        logLikelihood += Math.Log(sum);
      }
      return logLikelihood;
    }

    public static double computeBIC(double logLik) {
      return 2 * logLik - numClass * Math.Log(input.Length);
    }

    public static double normalDensity(double x, double mean, double stddev)
    {
      double diffSqr = Math.Pow((double)(x - mean), 2.0);
      double stddevSqr = 2 * Math.Pow((double)stddev, 2.0);
      double power = Math.Exp(-1.0 * diffSqr / stddevSqr);
      return (normalDensityConst * power);
    }

    public static double computeMean(int s, int e) {
      double sum = 0.0;
      for (int i = s; i < e; i++) {
        sum += input[i];
      }
      return sum/(e - s);
    }

    public static double computeStddev(int s, int e, int meanIdx) {
      double stddev = 0.0;
      for (int i = s; i < e; i++) {
        stddev += Math.Pow(input[i] - meanArr[0][meanIdx],2);
      }
      stddev/= (e - s - 1);
      stddev = Math.Sqrt(stddev);
      return stddev;
    }

    private static void print(double[][] arr) {
      for (int i = 0; i < arr.Length; i++) {
        for (int j = 0; j < arr[0].Length; j++) {
          Console.Write(addTrailingSpace(arr[i][j]));
        }
        Console.WriteLine();
      }
    }

    private static string addTrailingSpace(double val, int maxLen = 35) {
      StringBuilder x = new StringBuilder(val.ToString(), maxLen);
      for (int i = 1; i <= (maxLen - x.Length); i++) {
        x.Append(" ");
      }

      return x.ToString();
    }
  }
}
