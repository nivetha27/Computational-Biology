using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Assignment5
{
  class Program
  {
    public static int countASeqEnd = 10;
    static void Main(string[] args)
    {
      string A = null;
      for (int i = 1; i <= countASeqEnd; i++) {
        A += 'A';
      }

      using (StreamReader sr = new StreamReader(@"C:\Users\nsathya\Documents\GitHub\Computational-Biology\Assignment5\Data\SRR5831944.resorted2.sam")) {
        int totalRows = 0;
        int validRecords = 0;
        string s;
        while(!String.IsNullOrEmpty(s = sr.ReadLine())) {
          s = s.Trim();
          if (s.StartsWith("@")) // beginning of sam file
            continue;
        
          totalRows += 1;
          var sam = new Sam(s);
          if (sam.AS == null || sam.NM == null || sam.MDZ == null) {
            continue;
          }

          if (sam.ssegmentSeq.EndsWith(A)) {
            int putativeClevageIndex;
            var percentage = PercentageMismatchesInPolyTail(sam, out putativeClevageIndex);
            int polyTailLength = sam.ssegmentSeq.Length - putativeClevageIndex;
            //if (polyTailLength > 60)
            //{
            //    continue;
            //}

            // if (sam.NM > polyTailLength * 0.7 && percentage > 0.7 && sam.AS < -10)
            if (sam.AS <= -3 && sam.NM >= 5)
            {
              validRecords += 1;

              Console.WriteLine("[{0}]\t{1}\t{2}\t{3}\t{4}\tAS:{5}\tNM:{6}", validRecords, sam.queryTemplateName, sam.refSeqname, sam.pos, sam.cigar, sam.AS, sam.NM);
              Console.WriteLine("MD: {0}", sam.MDZ);
              Console.WriteLine("Putative cleavage location: {0}", putativeClevageIndex + 1);
              Console.WriteLine("SEQ: {0} {1}", sam.ssegmentSeq.Substring(0, putativeClevageIndex + 1), sam.ssegmentSeq.Substring(putativeClevageIndex + 1));
              Console.WriteLine();
            }
          }
        }
        Console.Write("{0} Valid Records out of {1}", validRecords, totalRows);
      }
    }
  
    
    static double PercentageMismatchesInPolyTail(Sam record, out int putativeClevageIndex)
    {
        putativeClevageIndex = record.ssegmentSeq.Length - 1;
        var mdCode = record.MDZ;

        int matchCount = 0;
        int mismatch = 0;
        int endGenome = EndGenome(record.ssegmentSeq);
        int units = 0;

        for (int i = mdCode.Length - 1; i >= 0 && putativeClevageIndex >= endGenome; i--)
        {
            var current = mdCode[i];
            if (char.IsDigit(current))
            {
                matchCount = int.Parse(current.ToString()) * (int)Math.Pow(10, units++) + matchCount;
            }
            else if (char.IsLetter(current) && (putativeClevageIndex - matchCount) > endGenome)
            {
                putativeClevageIndex -= matchCount;

                matchCount = 0;
                units = 0;

                if (putativeClevageIndex > endGenome)
                {
                    mismatch++;
                    putativeClevageIndex--;
                }
            }
        }

        return (double)mismatch / (int)record.NM;
    }

    static int EndGenome(string sequence)
    {
        int index = sequence.Length - 1;
        while (sequence[index] == 'A')
        {
            index--;
        }

        return index;
    }  
  }
}
