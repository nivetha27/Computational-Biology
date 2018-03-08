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
    public const int motifLen = 6;
    public const int numNucleotides = 4;
    public const string seqLabel = "SEQ: ";
    public const char cleavageIndicator = '.';
    public enum NucleotideEnum {A,C,G,T};

    public static string consensus = "AATAAA";

    public static int countASeqEnd = 6;
    public static int alignmentScore = -1;
    public static string inputSamFile = @"C:\Users\nsathya\Documents\GitHub\Computational-Biology\Assignment5\Data\SRR5831944.resorted2.sam";
    public static string candidateFile = @".\..\..\candidates.txt";
    
    public static List<string> candidates = new List<string>();

    static void Main(string[] args)
    {
      bool runFindCandidates = false;
      bool runToyExample = false;
      if (runToyExample) {
        candidates.Add("CTTCGAAGCGAAAAGTCCTAATAGTAGAAGAACCCTCCATAAACCTGGAGTGACTATATGGATGCCCCCCACCCTACCACACATTCGAAGAAC");
        candidates.Add("CTCAAAAAAAAAAAAAAAAGATAATGGCTTCTTGAAAAAACAAAGAAATCAACCTGAAGGAATTCCTGATGGCCAAAGCTAGAACAATCTGAG");
        candidates.Add("CGGTTTAAGAATACATCCTTGTATAATCTGACATACAAATTTGTCATTTCCTGCACATGCACACCATTGTTAAAAAAAAAAAAAAAAAGCCAG");
      } else {
        if (!runFindCandidates && File.Exists(candidateFile)) {
          using (StreamReader sr = new StreamReader(candidateFile)) {
            string s;
            while (!sr.EndOfStream) {
              s = sr.ReadLine();
              if (s.StartsWith(seqLabel))
              {
              string line = s.Replace(seqLabel, null);
              candidates.Add(line.Substring(0, line.IndexOf(cleavageIndicator)));
              }
            }
          }
        } else {
          findCandidates();
        }
      }

      WMM();

      Console.WriteLine("Press any key to exit...");
      Console.ReadLine();
    }

    public static void WMM() {
      double[][] p0 = initializeFrequencies(1, 0);
      double[][] p1 = initializeFrequencies(0.85, 0.05);
      calculate(p0, "WMM0");
      calculate(p1, "WMM1");
      double[][] p2 = computeP2(candidates, p1);
      calculate(p2, "WMM2");
    }

    public static void findCandidates() {
      string A = new String('A', countASeqEnd); //null;

      using (StreamReader sr = new StreamReader(inputSamFile)) {
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
            if (sam.AS <= alignmentScore && sam.softClip >= countASeqEnd && sam.NM >= 0)
            {
              validRecords += 1;
              int putativeClevageIndex = sam.ssegmentSeq.Length - (int)sam.softClip;
              Console.WriteLine("[{0}]\t{1}\t{2}\t{3}\t{4}\tAS:{5}\tNM:{6}\tMD:{7}", validRecords, sam.queryTemplateName, sam.refSeqname, sam.pos, sam.cigar, sam.AS, sam.NM, sam.MDZ);
              Console.WriteLine("Putative cleavage location: {0}", putativeClevageIndex);
              string seq = String.Format("{2}{0}{3}{1}", sam.ssegmentSeq.Substring(0, putativeClevageIndex), sam.ssegmentSeq.Substring(putativeClevageIndex), seqLabel, cleavageIndicator);
              Console.WriteLine(seq);
              Console.WriteLine();
              candidates.Add(sam.ssegmentSeq.Substring(0, putativeClevageIndex));
            }
          }
        }
        Console.Write("{0} Valid Records out of {1}", validRecords, totalRows);
      }
    }
    
    public static double getRelativeEntropy(double[][] wmm, double[][] p) {
      double relativeEntropy = 0;
      for (int i =  0; i < wmm.Length; i++) {
        for (int j = 0; j < wmm[0].Length; j++) {
          if (p[i][j] != 0) {
            relativeEntropy += (p[i][j] * wmm[i][j]);
          }
        }
      }

      return relativeEntropy;
    }

    public static void calculate(double[][] p, string label) {
      Console.WriteLine(label);
      Console.WriteLine();

      var wmm = computeWMMFromFreq(p);

      Console.WriteLine("Probabilities");
      print(p);
      Console.WriteLine();

      Console.WriteLine(label);
      print(wmm);
      Console.WriteLine();

      int positiveLLR = 0;
      int cumulativeDistance = 0;

      Console.WriteLine("\tMotif\tScore\tLength");

      foreach (var read in candidates) {
        // compute best hit
        double bestHitScore = Double.MinValue;
        string bestHitMotif = null;
        for (int i = 0; i <= read.Length - motifLen; i++) {
          var motif = read.Substring(i, motifLen);
          double score = 0;
          for (int j = 0; j < motif.Length; j++) {
            score += wmm[getNucleotideIdx(motif[j])][j];
          }

          if (score >= bestHitScore) {
            bestHitScore = score;
            bestHitMotif = motif;
          }
        }


        if (bestHitScore > 0) {
          var distance = read.Length - read.LastIndexOf(bestHitMotif);
          cumulativeDistance += distance;

          Console.WriteLine("[{0}]\t{1}\t{2}\t{3}", ++positiveLLR, bestHitMotif, bestHitScore.ToString("G4"), distance);
        }
      }

      if (positiveLLR > 0) {
        Console.WriteLine("Candidate count = {0}", candidates.Count);
        Console.WriteLine("Candidate count with positive LLR = {0}", positiveLLR);
        Console.WriteLine("Average distance = {0}", (double)cumulativeDistance / positiveLLR);
      } else {
        Console.WriteLine("No match");
      }
      Console.WriteLine("Relative entropy = {0}", getRelativeEntropy(wmm, p));

      Console.WriteLine();
    }

    public static double[][] computeP2(List<string> reads, double[][] p) {
      int rowCount = p.Length;
      int colCount = p[0].Length;
      var aligned = new List<SixMer>();
      double totalWeight = 0;
      foreach (var read in reads) {
        var sixMers = new List<SixMer>();
        double sequenceTotalWeight = 0;
        for (int i = 0; i <= read.Length - motifLen; i++) {
          var sixMer = read.Substring(i, motifLen);
          var probability = getProbability(sixMer, p);
          sixMers.Add(new SixMer(sixMer, probability));
          sequenceTotalWeight += probability;
        }
        
        if (sequenceTotalWeight > 0) {
          foreach (var sm in sixMers) {
            aligned.Add(new SixMer(sm.seq, sm.weight / sequenceTotalWeight));
            totalWeight += sm.weight / sequenceTotalWeight;
          }
        }
      }

      var weightedFrequencies = new double[rowCount][];
      for (int i =  0; i < weightedFrequencies.Length; i++) {
        weightedFrequencies[i] = new double[colCount];
      }

      foreach (var sixMer in aligned) {
        for (int i = 0; i < sixMer.seq.Length; i++) {
          weightedFrequencies[getNucleotideIdx(sixMer.seq[i])][i] += sixMer.weight;
        }
      }

      for (int i = 0; i < rowCount; i++)
      {
        for (int j = 0; j < colCount; j++)
        {
          weightedFrequencies[i][j] /= totalWeight;
        }
      }

      return weightedFrequencies;
    }

    public static void print(double[][] arr) {
      for (int i = 0; i < arr.Length; i++)
      {
        for (int j = 0; j < arr[0].Length; j++)
        {
          Console.Write("{0}\t", arr[i][j].ToString("G4"));
        }
        Console.WriteLine();
      }
    }

    public static double getProbability(string sixMer, double[][] p)
    {
      double probability = 1;
      for (int i = 0; i < sixMer.Length; i++)
      {
        probability *= (p[getNucleotideIdx(sixMer[i])][i]);
      }

      return probability;
    }

    public static int getNucleotideIdx(char a) {
      NucleotideEnum aEnum;
      int x = -1;
      if(Enum.TryParse<NucleotideEnum>(a.ToString(), true, out aEnum))
      {
        x = (int)aEnum;
      }
      return x;
    }

    public static double[][] computeWMMFromFreq(double[][] freqArray) {
      int row = freqArray.Length;
      int col = freqArray[0].Length;
      var wmm = new double[row][];
      for (int i = 0; i < row; i++)
      {
        wmm[i] = new double[col];
        for (int j = 0; j < col; j++)
        {
          wmm[i][j] = Math.Log(freqArray[i][j] / 0.25, 2);
        }
      }

      return wmm;
    }  

    public static double[][] initializeFrequencies(double consensusFreq, double nonConsensusFreq) {
      var frequencies = new double [numNucleotides][];
      for (int i = 0; i < numNucleotides; i++) {
        frequencies[i] = new double[consensus.Length];
        for (int j = 0; j < consensus.Length; j++) {
          frequencies[i][j] = nonConsensusFreq;
          if (i == getNucleotideIdx(consensus[j])) {
            frequencies[i][j] = consensusFreq;
          }
        }
      }
      return frequencies;
    }
  }
}
