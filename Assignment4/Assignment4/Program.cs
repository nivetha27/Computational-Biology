﻿﻿using System;
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

  public class Coordinates {
    public int start;
    public int end;
    public Coordinates(int s = -1, int e = -1) {
      start = s;
      end = e;
    } 
  }
  class Program
  {
    public const int numCharsToPrintPerLine = 60;
    public const int numStates = 2;
    public const int state0 = 0;
    public const int state1 = 1;
    public const char defaultCharReplaceNonACGT = 'T';
    public static string fileDirectory = @"C:\Users\nsathya\Documents\GitHub\Computational-Biology\Assignment4\Assignment4\";
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
      sequences = readAFastaFile(fileDirectory + "GCF_000091665.1_ASM9166v1_genomic.fna");
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
            // rolls[row][col].value = (HMMParams[numStates, row]) * (emissionStates[c][row]);
            rolls[row][col].value = Math.Log(HMMParams[numStates, row]) + Math.Log(emissionStates[c][row]);
          } else {
            rolls[row][col].value = Double.MinValue;
            for (int i = 0; i < numStates; i++) {
              // double stateValues = (rolls[i][col - 1].value) * (HMMParams[i, row]) * (emissionStates[c][row]);
              double stateValues = (rolls[i][col - 1].value) + Math.Log(HMMParams[i, row]) + Math.Log(emissionStates[c][row]);
              rolls[row][col].value = Math.Max(rolls[row][col].value, stateValues);
              rolls[row][col].allStates[i] = stateValues;
            }
          }
          maxLenInRolls = Math.Max(maxLenInRolls, rolls[row][col].ToString().Length);
        }
      }
    }

    /// <summary>
    /// This function traces back and determines
    /// 1. Viterbi path and prints it
    /// 2. #Hits (sequences of 1) and prints the location and length.
    /// </summary>
    public static void traceBack(int numHitsToPrnt = Int32.MaxValue, bool shouldPrintViterbi = true) {
      bool curState1SeqFound = false;
      StringBuilder mostProbablePath = new StringBuilder();
      List<Coordinates> hits = new List<Coordinates>();
      Coordinates x = new Coordinates();
      int prevState = -1;
      for (int j = sequences.Length - 1; j >= 0; j--) {
        if (j == sequences.Length - 1) {
          prevState = (rolls[state0][j].value > rolls[state1][j].value)  ? state0 : state1;
          Console.WriteLine(String.Format("Overall log probability {0}", Math.Exp(rolls[prevState][j].value)));
        } else {
          int k = j + 1;
          char c = validateSequenceChar(sequences[k]);
          prevState = ((prevState == state1 && rolls[state1][k].value == rolls[state1][k].allStates[state1]) || 
                       (prevState == state0 && rolls[state0][k].value == rolls[state0][k].allStates[state1])) ? state1  : state0;
        }
        mostProbablePath.Append(prevState.ToString());
        // mostProbablePath.Append(prevState == state0 ? "L" : "F");

        // calculate hits
        if (prevState == state1) {
          if (!curState1SeqFound) {
            x = new Coordinates(j + 1,j + 1);
            curState1SeqFound = true;
          } else {
            x.start = j + 1;
          }
          if (j == 0) { // edge case if the first element is state1
            hits.Add(x);
          }
        } else if (prevState == state0) {
          if (curState1SeqFound) {
            hits.Add(x);
          }
          curState1SeqFound = false;
        }
      }

      printHits(hits, numHitsToPrnt);
      if (shouldPrintViterbi)
        printViterbiPath(reverse(mostProbablePath.ToString()));
    }

    public static void printViterbiPath(string viterbi) {
      // Console.WriteLine(String.Format("Viterbi path is: {0}", viterbi));
      using (StreamWriter sw = new StreamWriter(fileDirectory + "output" + DateTime.Now.Ticks + ".txt")) {
        for (int i = 0; i < sequences.Length;)
        {
          string prefixFirst = "Rolls";
          string prefixSecond = "Viterbi";
          prefixFirst = addTrailiingWhiteSpaces(prefixFirst, 9 - prefixFirst.Length);
          prefixSecond = addTrailiingWhiteSpaces(prefixSecond, 9 - prefixSecond.Length);

          //Console.WriteLine(prefixFirst + " " + safeSubstring(sequences, i, numCharsToPrintPerLine));
          //Console.WriteLine(prefixSecond + " " + safeSubstring(viterbi, i, numCharsToPrintPerLine));
          //Console.WriteLine();

          sw.WriteLine(prefixFirst + " " + safeSubstring(sequences, i, numCharsToPrintPerLine));
          sw.WriteLine(prefixSecond + " " + safeSubstring(viterbi, i, numCharsToPrintPerLine));
          sw.WriteLine();
          i += numCharsToPrintPerLine;
        }
      }
    }

    public static void printHits(List<Coordinates> hits, int numHitsToPrnt = Int32.MaxValue) {
      Console.WriteLine(String.Format("#hits : {0}", hits.Count()));

      if (numHitsToPrnt > hits.Count()) {
        numHitsToPrnt = hits.Count();
      }
      Console.WriteLine(String.Format("Printing {0} / {1} hits", numHitsToPrnt, hits.Count()));
      for (int i = 0; i < numHitsToPrnt; i++) {
        Console.WriteLine(String.Format("Location = {0}, Length = {1}", hits[i].start, hits[i].end - hits[i].start + 1));
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
         { 0.01,0.99}, // State 2
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
