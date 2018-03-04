using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace WMM
{
  public class Program
  {
    public const int numNucleotides = 4;
    public const string seqLabel = "SEQ: ";
    public const char cleavageIndicator = ' ';
    public enum NucleotideEnum {A,C,G,T};

    public static string fullInputFileName = @"C:\Users\nsathya\Documents\GitHub\Computational-Biology\Assignment5\WMM\candidates.txt";
    public static string consensus = "AATAAA";
        
    public static double[,] p0;
    public static double[,] p1;

    public static void Main(string[] args)
    {
      p0 = initializeFrequencies(1, 0);
      p1 = initializeFrequencies(0.85, 0.05);

      List<string> reads = new List<string>();
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
      calculate(reads, p0, "WMM0");
      calculate(reads, p1, "WMM1");
      calculate(reads, getP2(reads, p1), "WMM2");
    }

    public static double[,] initializeFrequencies(double consensusFreq, double nonConsensusFreq) {
      var frequencies = new double [numNucleotides, consensus.Length];
      for (int i = 0; i < numNucleotides; i++) {
        for (int j = 0; j < consensus.Length; j++) {
          frequencies[i, j] = nonConsensusFreq;
          if (i == getNucleotideIdx(consensus[j])) {
            frequencies[i, j] = consensusFreq;
          }
        }
      }
      return frequencies;
    }

    public static double[,] computeWMMFromFreq(double[,] freqArray)
    {
      var wmm = new double[freqArray.GetLength(0), freqArray.GetLength(1)];
      for (int i = 0; i < wmm.GetLength(0); i++)
      {
        for (int j = 0; j < wmm.GetLength(1); j++)
        {
          wmm[i, j] = Math.Log(freqArray[i, j] / 0.25, 2);
        }
      }

      return wmm;
    }

    public static void calculate(List<string> reads, double[,] p, string label)
    {
      Console.WriteLine(label);
      Console.WriteLine();

      var wmm = computeWMMFromFreq(p);

      Console.WriteLine("Probabilities");
      print(p);
      Console.WriteLine();

      Console.WriteLine("WMM");
      print(wmm);
      Console.WriteLine();

      int positiveLLR = 0;
      int cumulativeDistance = 0;

      Console.WriteLine("\tMotif\tScore\tLength");

      foreach (var read in reads)
      {
        // compute best hit
        double bestHitScore = Double.MinValue;
        string bestHitMotif = null;
        for (int i = 0; i < read.Length - 6; i++)
        {
          var motif = read.Substring(i, 6);
          double score = 0;
          for (int j = 0; j < motif.Length; j++)
          {
            score += wmm[getNucleotideIdx(motif[j]), j];
          }

          if (score >= bestHitScore)
          {
            bestHitScore = score;
            bestHitMotif = motif;
          }
        }


        if (bestHitScore > 0)
        {
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

    public static double getRelativeEntropy(double[,] wmm, double[,] p)
    {
      double relativeEntropy = 0;
      for (int j = 0; j < wmm.GetLength(1); j++)
      {
        for (int i = 0; i < wmm.GetLength(0); i++)
        {
          if (p[i, j] != 0)
          {
            relativeEntropy += (p[i, j] * wmm[i, j]);
          }
        }
      }

      return relativeEntropy;
    }

    public static double[,] getP2(List<string> reads, double[,] p)
    {
      var aligned = new List<Tuple<string, double>>();
      foreach (var read in reads)
      {
        var sixMers = new List<Tuple<string, double>>();
        for (int i = 0; i < read.Length - 6; i++)
        {
          var sixMer = read.Substring(i, 6);
          var probability = getProbability(sixMer, p);

          sixMers.Add(new Tuple<string, double>(sixMer, probability));
        }

        var sequenceTotalWeight = sixMers.Sum(s => s.Item2);
        if (sequenceTotalWeight > 0)
        {
          foreach (var sm in sixMers)
          {
            aligned.Add(new Tuple<string, double>(sm.Item1, sm.Item2 / sequenceTotalWeight));
          }
        }
      }

      var totalWeight = aligned.Sum(s => s.Item2);
      var weightedFrequencies = new double[p.GetLength(0), p.GetLength(1)];
      foreach (var sixMer in aligned)
      {
        for (int i = 0; i < sixMer.Item1.Length; i++)
        {
          weightedFrequencies[getNucleotideIdx(sixMer.Item1[i]), i] += sixMer.Item2;
        }
      }

      for (int i = 0; i < weightedFrequencies.GetLength(0); i++)
      {
        for (int j = 0; j < weightedFrequencies.GetLength(1); j++)
        {
          weightedFrequencies[i, j] /= totalWeight;
        }
      }

      return weightedFrequencies;
    }

    public static void print(double[,] arr)
    {
      for (int i = 0; i < arr.GetLength(0); i++)
      {
        for (int j = 0; j < arr.GetLength(1); j++)
        {
          Console.Write("{0}\t", arr[i, j].ToString("G4"));
        }
        Console.WriteLine();
      }
    }

    public static double getProbability(string sixMer, double[,] p)
    {
      double probability = 1;
      for (int i = 0; i < sixMer.Length; i++)
      {
        probability *= (p[getNucleotideIdx(sixMer[i]), i]);
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
  }
}
