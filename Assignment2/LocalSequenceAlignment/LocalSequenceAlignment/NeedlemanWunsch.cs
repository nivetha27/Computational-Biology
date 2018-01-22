namespace SequenceAlignment
{
  // This class implements global sequence alignment
  public class NeedlemanWunsch: SequenceAlignment
  {
    public NeedlemanWunsch(string rowStr, string columnStr): base(rowStr, columnStr) {}

    public override bool traceBackLoopCondition(Score x)
    {
      return x.rowIdx > 0 && x.colIdx > 0;
    }
    public override int sequenceAlignment(int gapPenalty) {   
      // initialize first col as zero
      for(int row = 0; row < m; row++) {
        score[row] = new Score[n];
        score[row][0] = new global::SequenceAlignment.Score(row * gapPenalty);
      }

      // initialize first row as zero
      for (int col = 0; col < n; col++) {
        score[0][col] = new global::SequenceAlignment.Score(col * gapPenalty);
      }
      
      for (int row = 1; row < m; row++) {
        for (int col = 1; col < n; col++) {
          int value = 0;
          int x = row - 1, y = col - 1;

          int blosumVal = Blosum62.getBlosum62Value(rowStr[row - 1], columnStr[col - 1]);
          value = score[row - 1][col - 1].value + blosumVal;

          if (value < score[row][col-1].value + gapPenalty) {
            value = score[row][col-1].value + gapPenalty;
            x = row;
            y = col -1;
          }
          if (value < score[row-1][col].value + gapPenalty) {
            value = score[row-1][col].value + gapPenalty;
            x = row - 1;
            y = col;
          }

          score[row][col] = new global::SequenceAlignment.Score(value);
          score[row][col].parent = score[x][y];
          score[row][col].rowIdx = row;
          score[row][col].colIdx = col;
        }
      }
      maxScore = score[m - 1][n - 1];
      return maxScore.value;
    }
  }
}
