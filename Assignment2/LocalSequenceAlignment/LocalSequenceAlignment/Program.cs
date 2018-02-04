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

      Console.WriteLine("Default directory path is {0}", FastaFileProcessor.directory);
      Console.WriteLine("Do you want to update this ? (y/n)");
      string s = Console.ReadLine();
      if (s.ToUpper() == "Y") {
        Console.WriteLine("Enter directory of input fasta files : ");
        FastaFileProcessor.directory = Console.ReadLine();
        Console.WriteLine("Directory path successfully updated....");
      }

      List<Protein> proteins = null;
      Console.WriteLine("Select an option \n 1) local sequence alignment \n 2) global sequence alignment \n 3) local and global \n 4) emperical p-value \n 5) none");
      int alignmentOption = Convert.ToInt32(Console.ReadLine());
      if (alignmentOption < 4) {
        Console.WriteLine("Assuming required files exists, by default will run the local sequence alignment using protein sequences - P15172	P17542	P10085	P16075	P13904	Q90477	Q8IU24	P22816	Q10574	O95363");
        Console.WriteLine("Enter All[A] or pair of identifiers(ex: P15172,P17542) to run for");
        s = Console.ReadLine();        
        if (s.Length > 3) {
          proteins = FastaFileProcessor.getAllProteinSeqFromFile(s.Split(new char [] {','}, StringSplitOptions.RemoveEmptyEntries));
          Console.WriteLine("Protein sequences fetched....");
        } else if (s.ToUpper() == "ALL" || s.ToUpper() == "A" || proteins == null) {
          proteins = FastaFileProcessor.getAllProteinSeqFromFile();
        }    

        if (alignmentOption == 1 || alignmentOption == 3) {
          Console.WriteLine("\nRunning local sequence Alignment ... ");
          Controller.runLocalSequence(proteins, false, 100000);
        }
        
        if (alignmentOption == 2 || alignmentOption == 3) {
          Console.WriteLine("\nRunning global sequence Alignment ... ");
          Controller.runGlobalSequence(proteins);
        }        
      }
      if (alignmentOption == 4) {
        Console.WriteLine("Enter identifiers as comma separated list with at least two identifiers (ex: P15172,P17542) ");
        s = Console.ReadLine();
        proteins = FastaFileProcessor.getAllProteinSeqFromFile(s.Split(new char [] {','}, StringSplitOptions.RemoveEmptyEntries));
        Console.WriteLine("Protein sequences fetched....");
        Console.WriteLine("Enter number of random trials needed: ");
        s = Console.ReadLine();
        Console.WriteLine("Computing p-value");
        Controller.computeEmpericalPValue(null, proteins[0].sequence, proteins[1].sequence, Convert.ToInt32(s));
      }
      
      Console.WriteLine("Press any key to exit....");
      Console.ReadLine();     
    }
  }
}
