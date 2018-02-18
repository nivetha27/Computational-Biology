using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Assignment4
{
  public class StateEstimation {
    public double value;
    public double[] allStates;

    public StateEstimation(int numStates) {
      value = 0;
      allStates = new double[numStates];
    }
  }
  class Program
  {
    public const int numCharsToPrintPerLine = 60;
    public static char defaultCharReplaceNonACGT = 'T';
    public const int numStates = 2;
    public const int state1 = 0;
    public const int state2 = 1;
    // public static string sequences = "315116246446644245311321631164152133625144543631656626566666651166453132651245636664631636663162326455236266666625151631222555441666566563564324364131513465146353411126414626253356366163666466232534413661661163252562462255265252266435353336233121625364414432335163243633665562466662632666612355245242";
    public static string sequences = "316664";
    public static double[,] HMMParams;
    public static Dictionary<char, double[]> emissionStates = new Dictionary<char, double[]>();
    public static int maxLenInRolls = Int32.MinValue;
    public static StateEstimation[][] rolls = new StateEstimation[numStates][];
    static void Main(string[] args)
    {
      // initializeTestTransitionAndEmissionStates();
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
        rolls[i] = new StateEstimation[sequences.Length];

      for (int col = 0; col < sequences.Length; col++) {
        for (int row = 0; row < numStates; row++) {
          char c = validateSequenceChar(sequences[col]);
          rolls[row][col] = new StateEstimation(numStates);
          if (col == 0) {
            rolls[row][col].value = (HMMParams[numStates, row]) * (emissionStates[c][row]);
            // rolls[row][col] = Math.Log(HMMParams[numStates, row]) + Math.Log(emissionStates[c][row]);
          } else {
            rolls[row][col].value = Double.MinValue;
            for (int i = 0; i < numStates; i++) {
              double stateValues = (rolls[i][col - 1].value) * (HMMParams[i, row]) * (emissionStates[c][row]);
              //double stateValues = Math.Log(rolls[i][col - 1]) + Math.Log(HMMParams[i, col]) + Math.Log(emissionStates[c][col]);
              rolls[row][col].value = Math.Max(rolls[row][col].value, stateValues);
              rolls[row][col].allStates[i] = stateValues;
            }
          }
          maxLenInRolls = Math.Max(maxLenInRolls, rolls[row][col].ToString().Length);
        }
      }
    }
    public static void traceBack() {
      StringBuilder mostProbablePath = new StringBuilder();
      int prevState = -1;
      for (int j = sequences.Length - 1; j >= 0; j--) {
        if (j == sequences.Length - 1) {
          prevState = (rolls[state1][j].value > rolls[state2][j].value)  ? state1 : state2; 
        } else {
          int k = j + 1;
          char c = validateSequenceChar(sequences[k]);
          prevState = ((prevState == state2 && rolls[state2][k].value == rolls[state2][k].allStates[state2]) || 
                       (prevState == state1 && rolls[state1][k].value == rolls[state1][k].allStates[state2])) ? state2  : state1;
        }
        mostProbablePath.Append(prevState.ToString());
        // mostProbablePath.Append(prevState == state1 ? "L" : "F");
      }
      string viterbi = reverse(mostProbablePath.ToString());
      // Console.WriteLine(String.Format("Viterbi path is: {0}", viterbi));
      for (int i = 0; i < sequences.Length;) {
        string prefixFirst = "Rolls";
        string prefixSecond = "Viterbi";
        prefixFirst = addTrailiingWhiteSpaces(prefixFirst, 9 - prefixFirst.Length);
        prefixSecond = addTrailiingWhiteSpaces(prefixSecond, 9 - prefixSecond.Length);

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
      //HMMParams = new double[numStates + 1, numStates] {
      //   {0.9, 0.10}, // {0.6, 0.40},
      //   {0.05, 0.95}, //  {0.17, 0.83},
      //   {0.52, 0.48} //begin
      //};

      HMMParams = new double[numStates + 1, numStates] {
         {0.6, 0.40},
         {0.17, 0.83},
         {0.52, 0.48} //begin
      };

      emissionStates.Add('1', new double[numStates] { 1.0 / 10, 1.0 / 6 }); // 0 = Loaded, 1 = Fair
      emissionStates.Add('2', new double[numStates] { 1.0 / 10, 1.0 / 6 });
      emissionStates.Add('3', new double[numStates] { 1.0 / 10, 1.0 / 6 });
      emissionStates.Add('4', new double[numStates] { 1.0 / 10, 1.0 / 6 });
      emissionStates.Add('5', new double[numStates] { 1.0 / 10, 1.0 / 6 });
      emissionStates.Add('6', new double[numStates] { 1.0 / 2, 1.0 / 6 });
    }

    public static void print(StateEstimation[][] a) {      
      for (int i = 0; i < a.Length; i++) {
        StringBuilder s = new StringBuilder(addTrailiingWhiteSpaces(i.ToString(), 1));
        for (int j =0; j < a[0].Length; j++) {
          string rollProbStr = a[i][j].value.ToString();
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
