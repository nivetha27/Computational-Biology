using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Assignment4
{
  class Program
  {
    public const int numStates = 2;
    public const int state1 = 0;
    public const int state2 = 1;
    public static string sequences = "316664";
    public static double[,] HMMParams = new double[numStates + 1, numStates] {
   //// {state 1, State 2}
   //   {.9999, .0001}, // State 1
   //   {0.01,0.99}, // State 2
   //   {.9999, .0001} // Begin  

   // test data 
   {0.60, 0.40},
   {0.17, 0.83},
   {0.52, 0.48}
    };
    public static Dictionary<char, double[]> emissionStates = new Dictionary<char, double[]>();
    public static int maxLenInRolls = Int32.MinValue;
    public static double[][] rolls = new double[numStates][];
    static void Main(string[] args)
    {
      initializeEmissionStates();
      Console.WriteLine(sequences);

      initializeTransitionEmissionArray();
      maxLenInRolls += 2;

      // trace back
      traceBack();
      print(rolls);

      Console.WriteLine("Press any key to exit.....");
      Console.ReadLine();
    }

    public static void initializeTransitionEmissionArray() {
      for (int i = 0; i < numStates; i++)
        rolls[i] = new double[sequences.Length];

      for (int j = 0; j < sequences.Length; j++) {
          if (j == 0) {
            rolls[state1][j] = HMMParams[numStates, state1] * emissionStates[sequences[j]][state1];
            rolls[state2][j] = HMMParams[numStates, state2] * emissionStates[sequences[j]][state2];
          } else {
            double state11 = rolls[state1][j-1] * HMMParams[state1,state1] * emissionStates[sequences[j]][state1];
            double state21 = rolls[state2][j-1] * HMMParams[state2,state1] * emissionStates[sequences[j]][state1];
            rolls[state1][j] = Math.Max(state11, state21);
            double state12 = rolls[state1][j-1] * HMMParams[state1,state2] * emissionStates[sequences[j]][state2];
            double state22 = rolls[state2][j-1] * HMMParams[state2,state2] * emissionStates[sequences[j]][state2];
            rolls[state2][j] = Math.Max(state12, state22);
          }

          maxLenInRolls = Math.Max(maxLenInRolls, Math.Max(rolls[state1][j].ToString().Length,rolls[state2][j].ToString().Length));
        }
    }

    public static void traceBack() {
      StringBuilder mostProbablePath = new StringBuilder();
      int prevState = -1;
      for (int j = sequences.Length - 1; j >= 0; j--) {
        if (j == sequences.Length - 1) {
          if (rolls[state1][j] > rolls[state2][j]) {
            mostProbablePath.Append(state1.ToString());
            prevState = state1;
          } else {
            mostProbablePath.Append(state2.ToString());
            prevState = state2;
          } 
        } else {
          int k = j + 1;
          double state21 = rolls[state2][k-1] * HMMParams[state2,state1] * emissionStates[sequences[k]][state1];
          double state22 = rolls[state2][k-1] * HMMParams[state2,state2] * emissionStates[sequences[k]][state2];
          if ((prevState == state2 && rolls[state2][k] == state22) || (prevState == state1 && rolls[state1][k] == state21)) {
            mostProbablePath.Append(state2.ToString());
            prevState = state2;
          } else {
            mostProbablePath.Append(state1.ToString());
            prevState = state1;
          }
        }
      }
      Console.WriteLine(String.Format("Viterbi path is: {0}", reverse(mostProbablePath.ToString())));
    }
    public static void initializeEmissionStates() {
      //emissionStates.Add('A', new double[numStates] {0.25, 0.20}); // {state 1, State 2}
      //emissionStates.Add('C', new double[numStates] {0.25, 0.30});
      //emissionStates.Add('G', new double[numStates] {0.25, 0.30});
      //emissionStates.Add('T', new double[numStates] {0.25, 0.20});

      emissionStates.Add('1', new double[numStates] {1.0/10, 1.0/6});
      emissionStates.Add('2', new double[numStates] {1.0/10, 1.0/6});
      emissionStates.Add('3', new double[numStates] {1.0/10, 1.0/6});
      emissionStates.Add('4', new double[numStates] {1.0/10, 1.0/6});
      emissionStates.Add('5', new double[numStates] {1.0/10, 1.0/6});
      emissionStates.Add('6', new double[numStates] {1.0/2, 1.0/6});
    }

    public static void print(double[][] a) {      
      for (int i = 0; i < a.Length; i++) {
        StringBuilder s = new StringBuilder(addTrailiingWhiteSpaces(i.ToString(), 1));
        for (int j =0; j < a[0].Length; j++) {
          string rollProbStr = a[i][j].ToString();
          s.Append(addTrailiingWhiteSpaces(rollProbStr, maxLenInRolls - rollProbStr.Length));
        }
        Console.WriteLine(s.ToString());
      }
    }

    public static string addTrailiingWhiteSpaces(string x, int c) {
      StringBuilder s = new StringBuilder(x);
      for (int i = 1; i <=c; i++)
        s.Append(" ");
       
       return s.ToString(); 
    }

    public static string reverse(string s) {
      char[] x = s.ToCharArray();
      Array.Reverse(x);
      return new string(x);
    } 
  }
}
