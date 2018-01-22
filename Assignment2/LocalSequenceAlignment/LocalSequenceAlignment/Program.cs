using System;
using System.Collections.Generic;

namespace SequenceAlignment
{
  class Program
  {
    
    static void Main(string[] args)
    {
      Console.WriteLine("\nRunning test data ... ");
      Controller.testData();

      List<Protein> proteins = FastaFileProcessor.getProteinSeqFromFile();
      Console.WriteLine("\nRunning local sequence Alignment ... ");
      Controller.runLocalSequence(proteins, false, 100000);
      Console.WriteLine("\nRunning global sequence Alignment ... ");
      Controller.runGlobalSequence(proteins);
    }
  }
}
