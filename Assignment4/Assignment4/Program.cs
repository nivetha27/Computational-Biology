using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Assignment4
{
  public class stateEstimation {
    public double value;
    public double[] allStates;

    public stateEstimation(int numStates) {
      value = 0;
      allStates = new double[numStates];
    }
  }
  class Program
  {
    public const int numCharsToPrintPerLine = 80;
    public static char defaultCharReplaceNonACGT = 'T';
    public const int numStates = 2;
    public const int state1 = 0;
    public const int state2 = 1;
    // public static string sequences = "315116246446644245311321631164152133625144543631656626566666651166453132651245636664631636663162326455236266666625151631222555441666566563564324364131513465146353411126414626253356366163666466232534413661661163252562462255265252266435353336233121625364414432335163243633665562466662632666612355245242";
    public static string sequences = "316664";
    public static double[,] HMMParams;
    public static Dictionary<char, double[]> emissionStates = new Dictionary<char, double[]>();
    public static int maxLenInRolls = Int32.MinValue;
    public static double[][] rolls = new double[numStates][];
    static void Main(string[] args)
    {
      initializeTransitionAndEmissionStates();
      sequences = readAFastaFile(@"C:\Users\nsathya\Documents\GitHub\Computational-Biology\Assignment4\Assignment4\GCF_000091665.1_ASM9166v1_genomic.fna");
      // Console.WriteLine(sequences);

      initializeTransitionEmissionArray();
      maxLenInRolls += 2;

      // trace back
      traceBack();
      // print(rolls);

      Console.WriteLine("Press any key to exit.....");
      Console.ReadLine();
    }

    public static char validateSequenceChar(char x){
      char c = x;
      if (!emissionStates.ContainsKey(c)) {
        c = defaultCharReplaceNonACGT;
      }
      return c;
    }

    public static void initializeTransitionEmissionArray() {
      for (int i = 0; i < numStates; i++)
        rolls[i] = new double[sequences.Length];

      for (int j = 0; j < sequences.Length; j++) {
        char c = validateSequenceChar(sequences[j]);
        if (j == 0) {
          rolls[state1][j] = Math.Log(HMMParams[numStates, state1]) + Math.Log(emissionStates[c][state1]);
          rolls[state2][j] = Math.Log(HMMParams[numStates, state2]) + Math.Log(emissionStates[c][state2]);
          //rolls[state1][j] = (HMMParams[numStates, state1]) * (emissionStates[c][state1]);
          //rolls[state2][j] = (HMMParams[numStates, state2]) * (emissionStates[c][state2]);
        } else {
          double state11 = Math.Log(rolls[state1][j - 1]) + Math.Log(HMMParams[state1, state1]) + Math.Log(emissionStates[c][state1]);
          double state21 = Math.Log(rolls[state2][j - 1]) + Math.Log(HMMParams[state2, state1]) + Math.Log(emissionStates[c][state1]);
          rolls[state1][j] = Math.Max(state11, state21);
          double state12 = Math.Log(rolls[state1][j - 1]) + Math.Log(HMMParams[state1, state2]) + Math.Log(emissionStates[c][state2]);
          double state22 = Math.Log(rolls[state2][j - 1]) + Math.Log(HMMParams[state2, state2]) + Math.Log(emissionStates[c][state2]);
          rolls[state2][j] = Math.Max(state12, state22);
          //double state11 = (rolls[state1][j - 1]) * (HMMParams[state1, state1]) * (emissionStates[c][state1]);
          //double state21 = (rolls[state2][j - 1]) * (HMMParams[state2, state1]) * (emissionStates[c][state1]);
          //rolls[state1][j] = Math.Max(state11, state21);
          //double state12 = (rolls[state1][j - 1]) * (HMMParams[state1, state2]) * (emissionStates[c][state2]);
          //double state22 = (rolls[state2][j - 1]) * (HMMParams[state2, state2]) * (emissionStates[c][state2]);
          //rolls[state2][j] = Math.Max(state12, state22);
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
            prevState = state1;
          } else {
            // mostProbablePath.Append("F"); 
            prevState = state2;
          } 
        } else {
          int k = j + 1;
          char c = validateSequenceChar(sequences[j]);
          double state21 = Math.Log(rolls[state2][k - 1]) + Math.Log(HMMParams[state2, state1]) + Math.Log(emissionStates[c][state1]);
          double state22 = Math.Log(rolls[state2][k - 1]) + Math.Log(HMMParams[state2, state2]) + Math.Log(emissionStates[c][state2]);
          //double state21 = (rolls[state2][k - 1]) * (HMMParams[state2, state1]) * (emissionStates[c][state1]);
          //double state22 = (rolls[state2][k - 1]) * (HMMParams[state2, state2]) * (emissionStates[c][state2]);
          if ((prevState == state2 && rolls[state2][k] == state22) || (prevState == state1 && rolls[state1][k] == state21)) {
            // mostProbablePath.Append("F");
            prevState = state2;
          } else {
            // mostProbablePath.Append("L");
            prevState = state1;
          }
        }
        mostProbablePath.Append(prevState.ToString());
      }
      string viterbi = reverse(mostProbablePath.ToString());
      Console.WriteLine(String.Format("Viterbi path is: {0}", viterbi));
      
      for (int i = 0; i < sequences.Length;) {
        string prefixFirst = String.Format("Rolls    ");
        string prefixSecond = String.Format("Viterbi  ");

        Console.WriteLine(prefixFirst + " " + safeSubstring(sequences, i, numCharsToPrintPerLine));
        Console.WriteLine(prefixSecond + " " + safeSubstring(viterbi, i, numCharsToPrintPerLine));
        Console.WriteLine();
        i += numCharsToPrintPerLine;
      }
    }

    private static string readAFastaFile(string fullFileName) {
      StringBuilder seq = new StringBuilder();
      bool firstSeqSeen = false;
      using (StreamReader sr = new StreamReader(fullFileName)) {
        String s = sr.ReadLine();
        while(s != null) {
          if(!String.IsNullOrEmpty(s.Trim())) {
            if (s.StartsWith(">")) {
              if (!firstSeqSeen) {
                firstSeqSeen = true;
              } else {
                break;
              }
            } else {
              seq.Append(s);
            }
          }
          s = sr.ReadLine();
        }
      }
      return seq.ToString();
    }
    public static void initializeTransitionAndEmissionStates() {
      HMMParams = new double[numStates + 1, numStates] {
      // {state 1, State 2}
         {.9999, .0001}, // State 1
         {0.01,0.99}, // State 2
         {.9999, .0001} // Begin 
      };
      emissionStates.Add('A', new double[numStates] { 0.25, 0.20 }); // {state 1, State 2}
      emissionStates.Add('C', new double[numStates] { 0.25, 0.30 });
      emissionStates.Add('G', new double[numStates] { 0.25, 0.30 });
      emissionStates.Add('T', new double[numStates] { 0.25, 0.20 });
    }

    public static void initializeTestTransitionAndEmissionStates() {


      HMMParams = new double[numStates + 1, numStates] {
         /*{0.9, 0.10},*/ {0.6, 0.40},
         /*{0.05, 0.95},*/  {0.17, 0.83},
         {0.52, 0.48} //begin
      };

      emissionStates.Add('1', new double[numStates] { 1.0 / 10, 1.0 / 6 }); // 0 = Loaded, 1 = Fair
      emissionStates.Add('2', new double[numStates] { 1.0 / 10, 1.0 / 6 });
      emissionStates.Add('3', new double[numStates] { 1.0 / 10, 1.0 / 6 });
      emissionStates.Add('4', new double[numStates] { 1.0 / 10, 1.0 / 6 });
      emissionStates.Add('5', new double[numStates] { 1.0 / 10, 1.0 / 6 });
      emissionStates.Add('6', new double[numStates] { 1.0 / 2, 1.0 / 6 });
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

    private static string safeSubstring(string text, int start, int length)
    {
        return text.Length <= start ? ""
            : text.Length - start <= length ? text.Substring(start)
            : text.Substring(start, length);
    }
  }
}
