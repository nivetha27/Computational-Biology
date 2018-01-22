namespace SequenceAlignment
{
  public class Score {
    public Score parent;
    public int value;
    public int rowIdx;
    public int colIdx;

    public Score(int value) {
      this.value = value; 
      this.parent = null;
      this.rowIdx = -1;
      this.colIdx = -1;
    }
  }
}
