using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Assignment3
{
  class Program
  {
    public static Random rand = new Random();
    public static double epsilon = 0.001;
    public static readonly int numIter = 11;
    static int numClass = 5;
    public static readonly double normalDensityConst = Math.Sqrt(1/(2 * Math.PI));
    public static string sample = "9 10 11 20 21 22 46 49 55 57"; // -6 -5 -4 0 4 5 6
    public static double[] input;
    public static double[] meanArr;
    public static double[] stddevArr;
    static double[][] normalDensityMatrix;
    static double[] nomralDensitySumMatrix;
    static double[][] expectedMatrix;
    static void Main(string[] args)
    { 
      input = sample.Split(' ').Select(n => double.Parse(n)).ToArray();
      for (int i = 1; i <= 5; i++) {
        Console.WriteLine(i);
        numClass = i;
        meanArr = new double[numClass];
        stddevArr = new double[numClass];
        initializeCluster();
        initialize();
        //for(int j = 0; j < expectedMatrix.Length; i++) {
        //  Console.Write(addTrailingSpace(input[j]));
        //  print(expectedMatrix[j]);
        //}
      }
      Console.WriteLine("Press any key to exit....");
      Console.ReadLine();
    }

    public static void initializeCluster() {
      var min = input.Min();
      var max = input.Max();

      for (int i = 0; i<numClass; i++)
      {
          meanArr[i] = rand.Next((int)min, (int)max);
          stddevArr[i] = 1.0;
      }
      //meanArr = new double[5] { 35,12,46,22,45}; 
      //stddevArr = new double[5] {1,1,1,1,1};
    }

    public static void initialize() {
      double? prevLikelihood = null;
      double? curLikelihood = null;
      do {
        if (curLikelihood != null) {
          prevLikelihood = Convert.ToDouble(curLikelihood);
        }
        // E-Step
        normalDensityMatrix = new double[input.Length][];
        nomralDensitySumMatrix = new double[input.Length];
        expectedMatrix = new double[input.Length][];

        for (int row = 0; row < input.Length; row++) {
          normalDensityMatrix[row] = new double[numClass];
          nomralDensitySumMatrix[row] = 0;
          for (int col = 0; col < numClass; col++) {
            normalDensityMatrix[row][col] =  normalDensity(input[row], meanArr[col],stddevArr[col]);
            nomralDensitySumMatrix[row] += normalDensityMatrix[row][col];
          }
        }
        for (int row = 0; row < input.Length; row++) {
          expectedMatrix[row] = new double[numClass];  
          for (int col = 0; col < numClass; col++) {
            expectedMatrix[row][col] =  normalDensityMatrix[row][col]/nomralDensitySumMatrix[row];
          }
        }
        print(meanArr);
        var l = computeLogLikelihood();
        curLikelihood = l;
        Console.Write(addTrailingSpace(l));
        Console.Write(computeBIC(Convert.ToDouble(curLikelihood)).ToString("N3"));
        Console.WriteLine();

        // M-Step
        meanArr = new double[numClass];
        for (int col = 0; col < numClass; col++) {
          double denom = 0;
          for (int row = 0; row < input.Length; row++) {
            meanArr[col] += expectedMatrix[row][col] * input[row];
            denom += expectedMatrix[row][col];
          }
          meanArr[col] /= denom;
        }
      } while (shouldTerminate(prevLikelihood, curLikelihood));
    }

    public static bool shouldTerminate(double? prev, double? cur) {
      if (prev == null || cur == null)
        return true;
      double diff = Convert.ToDouble(cur) - Convert.ToDouble(prev);
      return diff > epsilon ? true : false;
    }

    public static double computeLogLikelihood(double stddev = 1.0) {
      double logLikelihood = 0;
      for (int i = 0; i < normalDensityMatrix.Length; i++)
      {
        double sum = 0;
        for (int j = 0; j < normalDensityMatrix[0].Length; j++)
        {
          sum += (1.0 / numClass * normalDensityMatrix[i][j]);
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
        stddev += Math.Pow(input[i] - meanArr[meanIdx],2);
      }
      stddev/= (e - s - 1);
      stddev = Math.Sqrt(stddev);
      return stddev;
    }

    private static void print(double[] arr) {
      for (int j = 0; j < arr.Length; j++) {
        Console.Write(addTrailingSpace(arr[j]));
      }
      double logLikelihood = (computeLogLikelihood());
    }

    private static string addTrailingSpace(double val, int maxLen = 35) {
      StringBuilder x = new StringBuilder(val.ToString("N12"), maxLen);
      for (int i = 1; i <= (maxLen - x.Length); i++) {
        x.Append(" ");
      }

      return x.ToString();
    }
  }
}
