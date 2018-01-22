using System;
using System.Text;

namespace SequenceAlignment
{
  public abstract class SequenceAlignment
  {
  public string rowStr; // row
    public string columnStr; // col
    public int m;
    public int n;
    public Score[][] score;
    public Score maxScore;

    public SequenceAlignment(string rowStr, string columnStr) {
      this.rowStr = rowStr;
      this.columnStr = columnStr;
      m = rowStr.Length + 1;
      n = columnStr.Length + 1;
      score = new Score[m][];
      maxScore = new Score (-1);
    }

    public abstract int sequenceAlignment(int gapPenalty);

    public abstract bool traceBackLoopCondition(Score x);
    public void traceBack(string name1 = "N1", string name2 = "N2") {
      StringBuilder first = new StringBuilder();
      StringBuilder second = new StringBuilder();
      StringBuilder middle = new StringBuilder();
      //first.Append(rowStr[maxScore.rowIdx - 1]);
      //second.Append(columnStr[maxScore.colIdx - 1]);
      //middle.Append(this.blastChar(rowStr[maxScore.rowIdx - 1], columnStr[maxScore.colIdx - 1]));
      Score cur = maxScore; // maxScore.parent;
      int s1Idx = -1, s2Idx = -1;

      while (traceBackLoopCondition(cur)) {
        Score parent = cur.parent;
        s1Idx = cur.rowIdx;
        s2Idx = cur.colIdx;
        if (cur.rowIdx == parent.rowIdx && cur.colIdx != parent.colIdx) {
          first.Append("-");
          second.Append(columnStr[cur.colIdx - 1]);
          middle.Append(this.blastChar('-', columnStr[cur.colIdx - 1]));
        } else if (cur.colIdx == parent.colIdx && cur.rowIdx != parent.rowIdx) {
          first.Append(rowStr[cur.rowIdx - 1]);
          second.Append("-");
          middle.Append(this.blastChar(rowStr[cur.rowIdx - 1], '-'));
        } else if (cur.colIdx != parent.colIdx && cur.rowIdx != parent.rowIdx) {
          first.Append(rowStr[cur.rowIdx - 1]);
          second.Append(columnStr[cur.colIdx - 1]);
          middle.Append(this.blastChar(rowStr[cur.rowIdx - 1], columnStr[cur.colIdx - 1]));
        } else {
          throw new Exception("Invalid parent");
        }
        cur = cur.parent;
      }

      string firstStr = reverseString(first.ToString());
      string secondStr = reverseString(second.ToString());
      string middleStr = reverseString(middle.ToString());
      int numGapsFirst = 0;
      int numGapsSecond = 0;
      int maxPrefxLen = 10;
      for (int i = 0; i < firstStr.Length;) {
        string prefixFirst = String.Format("{0} {1}", name1, (s1Idx + i - numGapsFirst));
        string prefixSecond = String.Format("{0} {1}", name2 , (s2Idx + i - numGapsSecond));
        prefixFirst = addTrailingWhitespace(prefixFirst, maxPrefxLen - prefixFirst.Length);
        prefixSecond = addTrailingWhitespace(prefixSecond, maxPrefxLen - prefixSecond.Length);

        Console.WriteLine(prefixFirst + " " + safeSubstring(firstStr, i, 60));
        Console.WriteLine(addTrailingWhitespace("", maxPrefxLen) + " " + safeSubstring(middleStr, i, 60));
        Console.WriteLine(prefixSecond + " " + safeSubstring(secondStr, i, 60));
        Console.WriteLine();
        numGapsFirst += countGaps(safeSubstring(firstStr, i, 60));
        numGapsSecond += countGaps(safeSubstring(secondStr, i, 60));
        i += 60;
      }
    }

    /******************************HELPER FUNCTION*************************************************/
    public void printScore() {
      for (int i = 0; i < score.Length; i++) {
        for(int j = 0; j < score[0].Length; j++) {
          Console.Write(score[i][j].value + " ");
        }
        Console.WriteLine();
      }
    }

    public void printMaxValueDetails() {
      Console.WriteLine(maxScore.value + " " + rowStr[maxScore.rowIdx - 1] + " " + columnStr[maxScore.colIdx - 1]);
    }

    private int countGaps(string x, char gap = '-') {
        int count = 0;
        for (int i = 0; i < x.Length; i++) {
          if (x[i] == gap) {
            ++count;
          }
        }
        return count;
    }

    private string addTrailingWhitespace(string x, int count) {
      StringBuilder s = new StringBuilder(x);
      for(int i = 1; i <=count; i++) {
        s.Append(" ");
      }
      return s.ToString();
    }

    private string blastChar(char a, char b) {
      if (a == b) {
        return a.ToString();
      } else if (a == '-' || b == '-') {
        return " ";
      } else if (Blosum62.getBlosum62Value(a,b) > 0) {
        return "+";
      }
      return " ";
    }

    private string reverseString(string s) {
       char[] x = (s.ToString().ToCharArray());
       Array.Reverse(x);
       return new string(x);
    }

    private string safeSubstring(string text, int start, int length)
    {
        return text.Length <= start ? ""
            : text.Length - start <= length ? text.Substring(start)
            : text.Substring(start, length);
    }
  }
}
