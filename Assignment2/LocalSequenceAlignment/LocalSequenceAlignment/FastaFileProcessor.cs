using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace SequenceAlignment
{
  public static class FastaFileProcessor
  {
    public static string filepath = @"C:\Users\nsathya\Desktop\CSE527\Assignment2\v2\LocalSequenceAlignment\LocalSequenceAlignment\Data\input.txt";

    public static List<Protein> getProteinSeqFromFile() {
      List<Protein> proteins = new List<Protein>();
      string name = null; 
      StringBuilder seq = new StringBuilder();
      using (StreamReader sr = new StreamReader(FastaFileProcessor.filepath)) {
        String s = sr.ReadLine();
        while(s != null) {
          if(!String.IsNullOrEmpty(s.Trim())) {
            if (s.StartsWith(">")) {
              name = s.Split(new char[] {'|'})[1];
            } else {
              seq.Append(s);
            }
          } else {
            if (!String.IsNullOrEmpty(name)) {
              proteins.Add(new Protein(name, seq.ToString()));
            }
            name = null; 
            seq = new StringBuilder();
          }
          s = sr.ReadLine();
        }
      }
      return proteins;
    }
  }
}
