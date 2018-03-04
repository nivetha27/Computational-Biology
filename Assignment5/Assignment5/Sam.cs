using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Assignment5
{
  public class Sam
  {
    public string queryTemplateName;
    public int flag;
    public string refSeqname;
    public int pos; // 1-based leftmost mapping position
    public int mappingQuality;
    public string cigar;
    public string refNameOfNextRead;
    public int posOfNextRead;
    public int templateLen;
    public string ssegmentSeq;
    public string phredScaleQuality;
    public int? AS;
    public int? ZS;
    public int? XN;
    public int? XM;
    public int? XO;
    public int? XG;
    public int? NM;
    public string MDZ = null;
    public string YTZ = null;
    public int? YS;
    public int? NH;

    /*
     * input s should be a tab separated
     * */
    public Sam(string s, char delimiter = '\t') {
      var properties = s.Split(delimiter);

      queryTemplateName = properties[0];
      Int32.TryParse(properties[1], out flag);
      refSeqname = properties[2];
      Int32.TryParse(properties[3], out pos);
      Int32.TryParse(properties[4], out mappingQuality);
      cigar = properties[5];
      refNameOfNextRead = properties[6];
      Int32.TryParse(properties[7], out posOfNextRead);
      Int32.TryParse(properties[8], out templateLen);
      ssegmentSeq = properties[9];
      phredScaleQuality = properties[10];

      for (int i = 11; i < properties.Length; i++) {
        var p = properties[i];
        if (p.StartsWith("AS:I:", StringComparison.InvariantCultureIgnoreCase)) {
          AS = convertOrDefault(p, "AS:I:");
        } else if (p.StartsWith("ZS:I:", StringComparison.InvariantCultureIgnoreCase)) {
          ZS = convertOrDefault(p, "ZS:I:");
        } else if (p.StartsWith("XN:I:", StringComparison.InvariantCultureIgnoreCase)) {
          XN = convertOrDefault(p, "XN:I:");
        } else if (p.StartsWith("XM:I:", StringComparison.InvariantCultureIgnoreCase)) {
          XM = convertOrDefault(p, "XM:I:");
        } else if (p.StartsWith("XO:I:", StringComparison.InvariantCultureIgnoreCase)) {
          XO = convertOrDefault(p, "XO:I:");
        } else if (p.StartsWith("XG:I:", StringComparison.InvariantCultureIgnoreCase)) {
          XG = convertOrDefault(p, "XG:I:");
        } else if (p.StartsWith("NM:I:", StringComparison.InvariantCultureIgnoreCase)) {
          NM = convertOrDefault(p, "NM:I:");
        } else if (p.StartsWith("MD:Z:", StringComparison.InvariantCultureIgnoreCase)) {
          MDZ = p.ToUpper().Replace("MD:Z:", null);
        } else if (p.StartsWith("YT:Z:", StringComparison.InvariantCultureIgnoreCase)) {
          YTZ = p.ToUpper().Replace("YT:Z:", null);
        } else if (p.StartsWith("NH:I:", StringComparison.InvariantCultureIgnoreCase)) {
          NH = convertOrDefault(p, "NH:I:");
        } else if (p.StartsWith("YS:I:", StringComparison.InvariantCultureIgnoreCase)) {
          YS = convertOrDefault(p, "YS:I:");
        }  
      }
    }

    private int? convertOrDefault(string stringVal, string remove) {
      string x = stringVal;
      if (!String.IsNullOrEmpty(x))
        x = x.ToUpper().Replace(remove.ToUpper(), null);
      int tempVal;
      int? val = Int32.TryParse(x, out tempVal) ? tempVal : (int?)null;
      return val;
    }
  }
}