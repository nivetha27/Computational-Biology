﻿using System;
using System.Collections.Generic;
using System.Text;

namespace SequenceAlignment
{
  public static class Controller
  {
    public static Random rnd = new Random();
    public static int gapPenalty = -4; 
    public static int maxSeqSizeToPrintScore = 15;    

    public static void runLocalSequence(List<Protein> proteins, bool runpEmperical = false, int empericalDenom = 10) {
      int counter = 0; //expect count to be n * (n-1)/2
      int[][] similarity = new int[proteins.Count][];
      for (int i = 0; i < similarity.Length; i++)
      {
        similarity[i] = new int[proteins.Count];
        for (int j = 0; j < similarity.Length; j++)
        {
          similarity[i][j] = -1;
        }
      }

      for (int i = 0; i < proteins.Count - 1; i++)
      {
        for (int j = i + 1; j < proteins.Count; j++)
        {
          counter++;          
          SmithWaterman seq = new SmithWaterman(proteins[i].sequence, proteins[j].sequence);
          int optimalScore = seq.sequenceAlignment(gapPenalty);
          similarity[i][j] = optimalScore;
          Console.WriteLine(String.Format("{0}) Identifiers {1} and {2}", counter, proteins[i].name, proteins[j].name));
          Console.WriteLine(String.Format("Score {0}", optimalScore));
          Console.WriteLine(String.Format("Optimal Alignment"));
          seq.traceBack(proteins[i].name, proteins[j].name);
          if (proteins[i].sequence.Length < maxSeqSizeToPrintScore && proteins[j].sequence.Length < maxSeqSizeToPrintScore)
          {
            seq.printScore();
          }
          if (runpEmperical && proteins[i].name == "P15172" && (proteins[j].name == "Q10574" || proteins[j].name == "O95363"))
          {
            computeEmpericalPValue(optimalScore, proteins[i].sequence, proteins[j].sequence, empericalDenom);
          }
        }
      }

      showScore(proteins, similarity);
    }

     public static void runGlobalSequence(List<Protein> proteins, bool runpEmperical = false, int empericalDenom = 10) {
      int counter = 0; //expect count to be n * (n-1)/2
      int[][] similarity = new int[proteins.Count][];
      for (int i = 0; i < similarity.Length; i++)
      {
        similarity[i] = new int[proteins.Count];
        for (int j = 0; j < similarity.Length; j++)
        {
          similarity[i][j] = -1;
        }
      }

      for (int i = 0; i < proteins.Count - 1; i++)
      {
        for (int j = i + 1; j < proteins.Count; j++)
        {
          counter++;          
          NeedlemanWunsch seq = new NeedlemanWunsch(proteins[i].sequence, proteins[j].sequence);
          int optimalScore = seq.sequenceAlignment(gapPenalty);
          similarity[i][j] = optimalScore;
          Console.WriteLine(String.Format("{0}) Identifiers {1} and {2}", counter, proteins[i].name, proteins[j].name));
          Console.WriteLine(String.Format("Score {0}", optimalScore));
          Console.WriteLine(String.Format("Optimal Alignment"));
          seq.traceBack(proteins[i].name, proteins[j].name);
          if (proteins[i].sequence.Length < maxSeqSizeToPrintScore && proteins[j].sequence.Length < maxSeqSizeToPrintScore)
          {
            seq.printScore();
          }
          if (runpEmperical && proteins[i].name == "P15172" && (proteins[j].name == "Q10574" || proteins[j].name == "O95363"))
          {
            computeEmpericalPValue(optimalScore, proteins[i].sequence, proteins[j].sequence, empericalDenom);
          }
        }
      }

      showScore(proteins, similarity);
    }

    public static void testData(bool runLocalAlign = true) {
      string x1="deadly", x2="ddgearlyk";
      if (runLocalAlign) {
        SmithWaterman seq = new SmithWaterman(x1, x2);
        int optimalScore = seq.sequenceAlignment(gapPenalty);
        seq.printMaxValueDetails();
        seq.traceBack();
        seq.printScore();
        computeEmpericalPValue(optimalScore, x1, x2, 999);
      } else {
        NeedlemanWunsch seq = new NeedlemanWunsch(x1, x2);
        int optimalScore = seq.sequenceAlignment(gapPenalty);
        seq.printMaxValueDetails();
        seq.traceBack();
        seq.printScore();
        computeEmpericalPValue(optimalScore, x1, x2, 999);
      }
    }

    public static void computeEmpericalPValue(int? optimalScore, string input1, string input2, int N) {
      if(optimalScore == null) {
        SmithWaterman seq = new SmithWaterman(input1, input2);
        optimalScore = seq.sequenceAlignment(gapPenalty);
      }
      int k = 0;
      for (int i = 1 ; i <= N ; i++) {
        string permutedString = permuteString(input2);
        SmithWaterman seq = new SmithWaterman(input1, permutedString);
        int score = seq.sequenceAlignment(gapPenalty);
        if (score >= optimalScore) {
          ++k;
        }
      }

      Console.WriteLine(String.Format("emperical p-value for k {0}, N {1} = {2} ", k, N, (k + 1) * 1.0/(N + 1)));
    }
    
    private static void showScore(List<Protein> proteins, int[][] similarity) {
      const int nameSize = 6;
      Console.WriteLine("Showing score");
      Console.Write(addTrailingWhitespace(" ", 7) + proteins[0].name + "  ");
      for(int i = 1; i < proteins.Count; i++) {
        Console.Write(proteins[i].name + "  ");
      }
      Console.WriteLine();
      for (int i = 0; i < similarity.Length; i++)
      {
        Console.Write(proteins[i].name + "  ");
        for (int j = 0; j < similarity.Length; j++)
        {
          string tmp = similarity[i][j].ToString();
          if(j <= i) {
            tmp = " ";
          }         
          Console.Write(addTrailingWhitespace(tmp, nameSize - tmp.Length) + "  ");
        }
        Console.WriteLine();
      }
    }

    private static string permuteString(string s) {
      char[] c = s.ToCharArray();
      for(int i = c.Length -1 ; i> 0; i--) {
        int j = rnd.Next(0,i);
        char tmp = c[i];
        c[i] = c[j];
        c[j] = tmp;
      }

      return new string(c);
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
