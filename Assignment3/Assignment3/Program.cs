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
    public static string sample = "14.04 14.71 15.26 13.85 15.2 15.03 15.09 16.12 13.78 16.27 20.26 19.87 20.28 21.25 21.15 20.69 20.05 20.35 22.22 21.2 23.84 23.12 24.59 21.67 24.03 23.52 27.32 27.02 24.86 22.73";
    // public static string sample = "9 10 11 20 21 22 46 49 55 57"; // -6 -5 -4 0 4 5 6
    public static List<List<double>> meanLogLikBICList;
    public static double[] input;
    public static double[] meanArr;
    public static double[] stddevArr;
    static double[][] normalDensityMatrix;
    static double[] nomralDensitySumMatrix;
    static double[][] expectedMatrix;
    public static int maxMeanLength; // needed for printing;
    public static int maxExpectedLength; //needed for printing;
    public static double maxBicScore = Int32.MinValue;
    public static int kWithMaxBicScore;
    static void Main(string[] args)
    { 
      input = sample.Split(' ').Select(n => double.Parse(n)).ToArray();
      for (int i = 1; i <= 5; i++) {
        maxMeanLength = maxExpectedLength = Int32.MinValue;
        meanLogLikBICList = new List<List<double>>();

        Console.WriteLine(i);
        numClass = i;
        meanArr = new double[numClass];
        stddevArr = new double[numClass];
        // Initializes mean and stddev.
        initializeCluster();        

        // Runs EM.
        initialize();

        // Print data.
        maxMeanLength +=5;
        printMeanMatrix();
        maxExpectedLength +=5;
        printExpectedMatrix();
      }

      Console.WriteLine(String.Format("k = {0} has best BIC Score of {1}", kWithMaxBicScore, maxBicScore));
      Console.WriteLine("Press any key to exit....");
      Console.ReadLine();
    }

    public static void initializeCluster() {
      var min = input.Min();
      var max = input.Max();
      for (int i = 0; i < numClass; i++)
      {
        meanArr[i] = rand.Next((int)min, (int)max);
        stddevArr[i] = 1.0;
        maxMeanLength = Math.Max(maxMeanLength, meanArr[i].ToString().Length);
      }
      meanLogLikBICList.Add(meanArr.ToList<double>());
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
            maxExpectedLength = Math.Max(maxExpectedLength, expectedMatrix[row][col].ToString().Length);
          }
        }
        // print(meanArr, maxMeanLength);
        var l = computeLogLikelihood();
        curLikelihood = l;
        meanLogLikBICList.Last().Add(l);
        maxMeanLength = Math.Max(maxMeanLength, l.ToString().Length);
        var b = computeBIC(l); 
        meanLogLikBICList.Last().Add(b);
        maxMeanLength = Math.Max(maxMeanLength, b.ToString().Length);
        if (maxBicScore < b) {
          maxBicScore = b;
          kWithMaxBicScore = meanArr.Length;
        }

        // M-Step
        meanArr = new double[numClass];
        for (int col = 0; col < numClass; col++) {
          double denom = 0;
          for (int row = 0; row < input.Length; row++) {
            meanArr[col] += expectedMatrix[row][col] * input[row];
            denom += expectedMatrix[row][col];
          }
          meanArr[col] /= denom;
          maxMeanLength = Math.Max(maxMeanLength, meanArr[col].ToString().Length);
        }
        meanLogLikBICList.Add(meanArr.ToList<double>());
      } while (shouldTerminate(prevLikelihood, curLikelihood));
      meanLogLikBICList.Remove(meanLogLikBICList.Last());
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

    private static void printMeanMatrix() {
      for (int j = 0; j < meanArr.Length; j++) {
        var x = "mu" + (j+1);
        Console.Write(addTrailingWhitespace(x, maxMeanLength - x.Length));
      }
      Console.Write(addTrailingWhitespace("LogLik", maxMeanLength - 6));
      Console.Write(addTrailingWhitespace("BIC", maxMeanLength - 3));
      Console.WriteLine();
      for (int j = 0; j < meanLogLikBICList.Count(); j++)
      { 
        print(meanLogLikBICList[j].ToArray<double>(), maxMeanLength);
        Console.WriteLine();
      }
    }

    private static void printExpectedMatrix() {
      Console.Write(addTrailingWhitespace(" ", maxExpectedLength));
      Console.Write(addTrailingWhitespace("xi", maxExpectedLength - "xi".Length));
      for (int j =1; j <= numClass; j++) {
          var x = String.Format("P(cls {0} | xi)", (j).ToString());
          Console.Write(addTrailingWhitespace(x, maxExpectedLength - x.Length));
      }
      Console.WriteLine();

      for (int j = 0; j < expectedMatrix.Length; j++)
      { 
        var x = String.Format("[{0},]", (j + 1).ToString());
        Console.Write(addTrailingWhitespace(x, maxExpectedLength - x.Length));
        Console.Write(addTrailingWhitespace(input[j].ToString(), maxExpectedLength - input[j].ToString().Length));
        print(expectedMatrix[j], maxExpectedLength);
        Console.WriteLine();
      }
    }

    private static void print(double[] arr, int count) {
      for (int j = 0; j < arr.Length; j++) {
        Console.Write(addTrailingWhitespace(arr[j].ToString(), count - arr[j].ToString().Length));
      }
    }

    private static string addTrailingWhitespace(string x, int count) {
      StringBuilder s = new StringBuilder(x);
      for(int i = 1; i <=count; i++) {
        s.Append(" ");
      }
      return s.ToString();
    }
  }
}
