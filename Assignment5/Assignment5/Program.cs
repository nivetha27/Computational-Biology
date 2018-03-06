﻿using System;
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
    public const char cleavageIndicator = ' ';
    public enum NucleotideEnum {A,C,G,T};

    public static string consensus = "AATAAA";

    public static int countASeqEnd = 10;
    public static string inputSamFile = @"C:\Users\nsathya\Documents\GitHub\Computational-Biology\Assignment5\Data\SRR5831944.resorted2.sam";
    public static string fullInputFileName = @"C:\Users\nsathya\Documents\GitHub\Computational-Biology\Assignment5\WMM\candidates.txt";
    
    public static List<string> reads = new List<string>();

    static void Main(string[] args)
    {
      bool runFindCandidates = false;
      if (!runFindCandidates && File.Exists(fullInputFileName)) {
        using (StreamReader sr = new StreamReader(fullInputFileName)) {
          string s;
          while (!sr.EndOfStream) {
            s = sr.ReadLine();
            if (s.StartsWith(seqLabel))
            {
            string line = s.Replace(seqLabel, null);
            reads.Add(line.Substring(0, line.IndexOf(cleavageIndicator)));
            }
          }
        }
      } else {
        findCandidates();
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
      double[][] p2 = computeP2(reads, p1);
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
            int putativeClevageIndex;
            var percentage = percentageMismatchesInPolyTail(sam, out putativeClevageIndex);
            int polyTailLength = sam.ssegmentSeq.Length - putativeClevageIndex;
            //if (polyTailLength > 60)
            //{
            //    continue;
            //}

            // if (sam.NM > polyTailLength * 0.7 && percentage > 0.7 && sam.AS < -10)
            // if (sam.AS <= -3 && sam.NM >= 5)
            if (sam.softClip >= 1)
            {
              validRecords += 1;

              Console.WriteLine("[{0}]\t{1}\t{2}\t{3}\t{4}\tAS:{5}\tNM:{6}", validRecords, sam.queryTemplateName, sam.refSeqname, sam.pos, sam.cigar, sam.AS, sam.NM);
              Console.WriteLine("MD: {0}", sam.MDZ);
              Console.WriteLine("Putative cleavage location: {0}", putativeClevageIndex + 1);
              string seq = String.Format("{2}{0} {1}", sam.ssegmentSeq.Substring(0, putativeClevageIndex + 1), sam.ssegmentSeq.Substring(putativeClevageIndex + 1), seqLabel);
              Console.WriteLine(seq);
              Console.WriteLine();
              reads.Add(seq);
            }
          }
        }
        Console.Write("{0} Valid Records out of {1}", validRecords, totalRows);
      }
    }
    
    public static double percentageMismatchesInPolyTail(Sam record, out int putativeClevageIndex) {
        putativeClevageIndex = record.ssegmentSeq.Length - 1;
        var mdCode = record.MDZ;

        int matchCount = 0;
        int mismatch = 0;
        int endGenome = Program.endGenome(record.ssegmentSeq);
        int units = 0;
        bool isCharSeen = false;
        for (int i = mdCode.Length - 1; i >= 0 && putativeClevageIndex >= endGenome; i--)
        {
            var current = mdCode[i];
            if (char.IsDigit(current))
            {
                matchCount = int.Parse(current.ToString()) * (int)Math.Pow(10, units++) + matchCount;
            }
            else if (char.IsLetter(current) && (putativeClevageIndex - matchCount) > endGenome)
            {
                isCharSeen = true;
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


        //if (!isCharSeen) {
        //  putativeClevageIndex = matchCount;
        //  if (putativeClevageIndex > endGenome)
        //  {
        //      mismatch++;
        //      putativeClevageIndex--;
        //  }
        //}

        return (double)mismatch / (int)record.NM;
    }

    public static int endGenome(string sequence) {
        int index = sequence.Length - 1;
        while (index >= 0 && sequence[index] == 'A') {
            index--;
        }
        return index;
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

      foreach (var read in reads) {
        // compute best hit
        double bestHitScore = Double.MinValue;
        string bestHitMotif = null;
        for (int i = 0; i < read.Length - motifLen; i++) {
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

      Console.WriteLine("Candidate count = {0}", reads.Count);
      Console.WriteLine("Candidate count with positive LLR = {0}", positiveLLR);
      Console.WriteLine("Average distance = {0}", cumulativeDistance / positiveLLR);
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
        for (int i = 0; i < read.Length - motifLen; i++) {
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

      for (int i = 0; i < rowCount; i++) {
        for (int j = 0; j < colCount; j++) {
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
