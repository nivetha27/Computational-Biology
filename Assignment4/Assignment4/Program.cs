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
    public const double epsilon = 0.001;
    public const int maxiteration = 10;
    public const int numCharsToPrintPerLine = 60;
    public const int numStates = 2;
    public const int state0 = 0;
    public const int state1 = 1;
    public const char defaultCharReplaceNonACGT = 'T';
    public static string fileDirectory = @"C:\Users\nsathya\Documents\GitHub\Computational-Biology\Assignment4\Assignment4\";
    public static string outputFullFileName = fileDirectory + "output" + DateTime.Now.Ticks + ".txt";
    public static StreamWriter sw = new StreamWriter(outputFullFileName, append: true);
    public static string sequences = "316664";
    public static double[,] transitionStates;
    public static Dictionary<char, double[]> emissionStates = new Dictionary<char, double[]>();
    public static int maxLenInRolls = Int32.MinValue;
    public static StateEstimation[][] HMM = new StateEstimation[numStates][];
    static void Main(string[] args)
    {
      bool runTestData = true;
      if (runTestData) {
        initializeTestTransitionAndEmissionStates(toyExample: true);
      } else {
        initializeTransitionAndEmissionStates();
        sequences = readAFastaFile(fileDirectory + "GCF_000091665.1_ASM9166v1_genomic.fna");
        // Console.WriteLine(sequences);
      }
      // EMWithEpsilon();
      EMWithMaxIterations();      

      sw.Flush();
      sw.Close();
      Console.WriteLine("Press any key to exit.....");
      Console.ReadLine();
    }

    public static void EMWithMaxIterations() {
      for (int i = 1; i <= maxiteration; i++) {
        computeHMMForward();
        printTransitionAndEmissionProb();        
        string viterbi = traceBack();
        int viterbiLastIdx = viterbi.Length - 1;
        double curMaxLogProb = HMM[(int)Char.GetNumericValue(viterbi[viterbiLastIdx])][viterbiLastIdx].value;
        printToConsoleAndWriteToFile(String.Format("Overall log probability {0}", curMaxLogProb));
        int x = i < maxiteration ? maxiteration/2 : Int32.MaxValue;
        calculateHits(viterbi,x);
        printViterbiPath(viterbi);
        trainData(viterbi);
        // print(HMM);
      }
    }

    public static void EMWithEpsilon() {
      double curMaxLogProb = 0.0, prevMaxLogProb = 0.0;
      do {
        prevMaxLogProb = curMaxLogProb;
        computeHMMForward();
        printTransitionAndEmissionProb();
        
        string viterbi = traceBack();
        int viterbiLastIdx = viterbi.Length - 1;
        curMaxLogProb = HMM[(int)Char.GetNumericValue(viterbi[viterbiLastIdx])][viterbiLastIdx].value;
        printToConsoleAndWriteToFile(String.Format("Overall log probability {0}", curMaxLogProb));
        calculateHits(viterbi);
        printViterbiPath(viterbi);
        trainData(viterbi);
        // print(HMM);
      } while(Math.Abs(curMaxLogProb - prevMaxLogProb) > epsilon);
    }

    public static void computeHMMForward() {
      for (int i = 0; i < numStates; i++)
        HMM[i] = new StateEstimation[sequences.Length];

      for (int col = 0; col < sequences.Length; col++) {
        for (int row = 0; row < numStates; row++) {
          char c = validateSequenceChar(sequences[col]);
          HMM[row][col] = new StateEstimation(numStates);
          if (col == 0) {
            // HMM[row][col].value = (transitionStates[numStates, row]) * (emissionStates[c][row]);
            HMM[row][col].value = Math.Log(transitionStates[numStates, row]) + Math.Log(emissionStates[c][row]);
          } else {
            HMM[row][col].value = Double.MinValue;
            for (int i = 0; i < numStates; i++) {
              // double stateValues = (HMM[i][col - 1].value) * (transitionStates[i, row]) * (emissionStates[c][row]);
              double stateValues = (HMM[i][col - 1].value) + Math.Log(transitionStates[i, row]) + Math.Log(emissionStates[c][row]);
              HMM[row][col].value = Math.Max(HMM[row][col].value, stateValues);
              HMM[row][col].allStates[i] = stateValues;
            }
          }
          maxLenInRolls = Math.Max(maxLenInRolls, HMM[row][col].ToString().Length);
        }
      }
      maxLenInRolls += 2;
    }

    /// <summary>
    /// This function traces back and determines
    /// 1. Viterbi path and prints it
    /// 2. Prints the overall log probability
    /// </summary>
    public static string traceBack() {
      StringBuilder mostProbablePath = new StringBuilder();
      int prevState = -1;
      for (int j = sequences.Length - 1; j >= 0; j--) {
        if (j == sequences.Length - 1) {
          prevState = (HMM[state0][j].value > HMM[state1][j].value)  ? state0 : state1;
        } else {
          int k = j + 1;
          char c = validateSequenceChar(sequences[k]);
          prevState = ((prevState == state1 && HMM[state1][k].value == HMM[state1][k].allStates[state1]) || 
                       (prevState == state0 && HMM[state0][k].value == HMM[state0][k].allStates[state1])) ? state1  : state0;
        }
        mostProbablePath.Append(prevState.ToString());
        // mostProbablePath.Append(prevState == state0 ? "L" : "F");
      }

      return reverse(mostProbablePath.ToString());
    }

    public static void calculateHits(string viterbi, int k = Int32.MaxValue) {
      bool curState1SeqFound = false;
      List<Coordinates> hits = new List<Coordinates>();
      Coordinates x = new Coordinates();
      for (int j = 0; j < viterbi.Length ; j++) {
        int curState = (int)Char.GetNumericValue(viterbi[j]);
        if (curState == state1) {
          if (!curState1SeqFound) {
            x = new Coordinates(j + 1,j + 1);
            curState1SeqFound = true;
          } else {
            x.end = j + 1;
          }
          if (j == viterbi.Length - 1) { // edge case if the first element is state1
            hits.Add(x);
          }
        } else if (curState == state0) {
          if (curState1SeqFound) {
            hits.Add(x);
          }
          curState1SeqFound = false;
        }
      }
      printHits(hits, k);
    }

    public static void trainData(string viterbi) {
      double[,] newTransitionStates = new double[numStates + 1, numStates];
      newTransitionStates[numStates, state0] = transitionStates[numStates, state0];
      newTransitionStates[numStates, state1] = transitionStates[numStates, state1];

      Dictionary<char, double[]> newEmissionStates = new Dictionary<char, double[]>();
      double[] newEmissionStatesOverallCount = new double[numStates];   
      foreach (var emissionState in emissionStates) {
        newEmissionStates.Add(emissionState.Key, new double[numStates]);
      }

      // get all transition counts and emission counts.
      for (int i = 0; i < viterbi.Length; i++) {
        int curState = (int)Char.GetNumericValue(viterbi[i]);
        if (i > 0) {
          int prevState = (int)Char.GetNumericValue(viterbi[i-1]);
          newTransitionStates[prevState, curState] += 1;        
        }
        
        char x = validateSequenceChar(sequences[i]);
        if (newEmissionStates.ContainsKey(x)) {
          newEmissionStates[x][curState] += 1;
          newEmissionStatesOverallCount[curState] += 1;
        }
      }

      // dividing by the appropriate counts
      foreach(var states in newEmissionStates.Values) {
        for (int i = 0; i < states.Length; i++) {
          states[i] /= newEmissionStatesOverallCount[i];
        }
      }

      for (int i = 0; i < newTransitionStates.GetLength(0) - 1; i++) {
        double sum = 0.0;
        for (int j= 0; j < newTransitionStates.GetLength(1); j++) {
          sum += newTransitionStates[i, j];
        }
        for (int j= 0; j < newTransitionStates.GetLength(1); j++) {
          newTransitionStates[i, j] /= sum;
        }
      }

      transitionStates = newTransitionStates;
      emissionStates = newEmissionStates;
    }

    public static void initializeTransitionAndEmissionStates() {
      transitionStates = new double[numStates + 1, numStates] {
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

    public static void initializeTestTransitionAndEmissionStates(bool toyExample = true) {
      if (toyExample) {
        sequences = "316664";
        transitionStates = new double[numStates + 1, numStates] {
           {0.6, 0.40},
           {0.17, 0.83},
           {0.52, 0.48} //begin
        };
      } else {
        sequences = "315116246446644245311321631164152133625144543631656626566666651166453132651245636664631636663162326455236266666625151631222555441666566563564324364131513465146353411126414626253356366163666466232534413661661163252562462255265252266435353336233121625364414432335163243633665562466662632666612355245242";
        transitionStates = new double[numStates + 1, numStates] {
           {0.9, 0.10}, // {0.6, 0.40},
           {0.05, 0.95}, //  {0.17, 0.83},
           {0.52, 0.48} //begin
        };
      }
      emissionStates.Add('1', new double[numStates] { 1.0 / 10, 1.0 / 6 }); // 0 = Loaded, 1 = Fair
      emissionStates.Add('2', new double[numStates] { 1.0 / 10, 1.0 / 6 });
      emissionStates.Add('3', new double[numStates] { 1.0 / 10, 1.0 / 6 });
      emissionStates.Add('4', new double[numStates] { 1.0 / 10, 1.0 / 6 });
      emissionStates.Add('5', new double[numStates] { 1.0 / 10, 1.0 / 6 });
      emissionStates.Add('6', new double[numStates] { 1.0 / 2, 1.0 / 6 });
    }

    /*************************************************************************PRINT FUNCTIONS**********************************************************************************************************/
    public static void printToConsoleAndWriteToFile(string s, bool writeInNewline = true, bool writeToConsole = true) {
      if (String.IsNullOrEmpty(s)) {
        if (writeInNewline) {
          s = Environment.NewLine;
        }
      } else if (writeInNewline) {
        s = s + Environment.NewLine;
      }

      if (writeToConsole) {
        Console.Write(s);
      }
      sw.Write(s);
    }
    
    public static void printTransitionAndEmissionProb() {
      printToConsoleAndWriteToFile("Transition States Probabilities");
      printToConsoleAndWriteToFile(addTrailiingWhiteSpaces(" ", 10 - 1), false);
      for (int i = 0; i < numStates; i++) {
        string curState = "state" + (i+1).ToString();
        printToConsoleAndWriteToFile(addTrailiingWhiteSpaces(curState, maxLenInRolls - curState.Length), false);
      }
      printToConsoleAndWriteToFile(null);
      for (int i = 0; i < transitionStates.GetLength(0); i++) {
        string curState = i == transitionStates.GetLength(0) - 1 ? "begin" : "state" + (i+1).ToString();
        printToConsoleAndWriteToFile(addTrailiingWhiteSpaces(curState, 10 - curState.Length), false);
        for (int j = 0; j < transitionStates.GetLength(1); j++) {
          printToConsoleAndWriteToFile(addTrailiingWhiteSpaces(transitionStates[i,j].ToString(), maxLenInRolls - transitionStates[i,j].ToString().Length), false);
        }
        printToConsoleAndWriteToFile(null);
      }

      printToConsoleAndWriteToFile("Emission States Probabilities");
      printToConsoleAndWriteToFile(addTrailiingWhiteSpaces(" ", 4), false);
      for (int i = 0; i < numStates; i++) {
        string curState = "state" + (i+1).ToString();
        printToConsoleAndWriteToFile(addTrailiingWhiteSpaces(curState, maxLenInRolls - curState.Length), false);
      }
      printToConsoleAndWriteToFile(null);
      foreach(var emissionState in emissionStates) {
        printToConsoleAndWriteToFile(addTrailiingWhiteSpaces(emissionState.Key.ToString(), 4), false);
        for (int j = 0; j < emissionState.Value.Length; j++) {
          printToConsoleAndWriteToFile(addTrailiingWhiteSpaces(emissionState.Value[j].ToString(), maxLenInRolls - emissionState.Value[j].ToString().Length), false);
        }
        printToConsoleAndWriteToFile(null);
      }
    }
    
    public static void print(StateEstimation[][] a) {      
      for (int i = 0; i < a.Length; i++) {
        StringBuilder s = new StringBuilder(addTrailiingWhiteSpaces(i.ToString(), 1));
        for (int j =0; j < a[0].Length; j++) {
          string rollProbStr = a[i][j].value.ToString();
          s.Append(addTrailiingWhiteSpaces(rollProbStr, maxLenInRolls - rollProbStr.Length));
        }
        printToConsoleAndWriteToFile(s.ToString());
      }
    }

    public static void printViterbiPath(string viterbi) {
      for (int i = 0; i < sequences.Length;)
      {
        string prefixFirst = "Rolls";
        string prefixSecond = "Viterbi";
        prefixFirst = addTrailiingWhiteSpaces(prefixFirst, 9 - prefixFirst.Length);
        prefixSecond = addTrailiingWhiteSpaces(prefixSecond, 9 - prefixSecond.Length);

        printToConsoleAndWriteToFile(prefixFirst + " " + safeSubstring(sequences, i, numCharsToPrintPerLine), true, false);
        printToConsoleAndWriteToFile(prefixSecond + " " + safeSubstring(viterbi, i, numCharsToPrintPerLine), true, false);
        printToConsoleAndWriteToFile(null, true, false);
        i += numCharsToPrintPerLine;
      }
    }

    public static void printHits(List<Coordinates> hits, int numHitsToPrnt = Int32.MaxValue) {
      printToConsoleAndWriteToFile(String.Format("#hits : {0}", hits.Count()));

      if (numHitsToPrnt > hits.Count()) {
        numHitsToPrnt = hits.Count();
      }
      printToConsoleAndWriteToFile(String.Format("Printing {0} / {1} hits", numHitsToPrnt, hits.Count()));
      for (int i = 0; i < numHitsToPrnt; i++) {
        printToConsoleAndWriteToFile(String.Format("Start = {0}, End = {1}, Length = {2}", hits[i].start, hits[i].end, hits[i].end - hits[i].start + 1));
      }
    }

    /*************************************************************************HELPER FUNCTIONS********************************************************************************************************/
    private static char validateSequenceChar(char x){
      char c = x;
      if (!emissionStates.ContainsKey(c)) {
        c = defaultCharReplaceNonACGT;
      }
      return c;
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
    
    private static string addTrailiingWhiteSpaces(string x, int c) {
      StringBuilder s = new StringBuilder(x);
      for (int i = 1; i <=c; i++)
        s.Append(" ");
       
       return s.ToString(); 
    }

    private static string reverse(string s) {
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
