using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace SequenceAlignment
{
  public static class FastaFileProcessor
  {
    public static string directory = @"C:\Users\nsathya\Documents\GitHub\Computational-Biology\Assignment2\LocalSequenceAlignment\LocalSequenceAlignment\Data";
    public static string allProteinSeqFile = FastaFileProcessor.directory + "input.txt";
    public static string fastaFileFormat = ".fasta";

    public static List<Protein> getAllProteinSeqFromFile(string[] files = null) {
      if (files == null) {
        files = System.IO.Directory.GetFiles(FastaFileProcessor.directory, "*" + FastaFileProcessor.fastaFileFormat);
      } else {
        for (int i = 0; i < files.Length; i++) {
         files[i] = FastaFileProcessor.directory + @"\" + files[i].Trim() + FastaFileProcessor.fastaFileFormat;
        }
      }      
      List<Protein> proteins = new List<Protein>();
      foreach(string file in files) {
        proteins.Add(readAFastaFile(file));
      }
      return proteins;
    }

    private static Protein readAFastaFile(string fullFileName) {
      string name = null; 
      StringBuilder seq = new StringBuilder();
      using (StreamReader sr = new StreamReader(fullFileName)) {
        String s = sr.ReadLine();
        while(s != null) {
          if(!String.IsNullOrEmpty(s.Trim())) {
            if (s.StartsWith(">")) {
              name = s.Split(new char[] {'|'})[1];
            } else {
              seq.Append(s);
            }
          }
          s = sr.ReadLine();
        }
      }
      return new Protein(name, seq.ToString());
    }
  }
}
