using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NaiveShuffle
{
  class Program
  { 
    public static Random r = new Random();
    static void Main(string[] args)
    {
      Console.WriteLine("Naive");
      shuffle(1000000);
      Console.WriteLine("Fisher");
      fisher(1000000);
      Console.WriteLine("Exit");
      Console.ReadLine();
    }
    
    static void shuffle(int N) {      
      Dictionary<string, int> map = new Dictionary<string, int>();

      for (int i = 0; i < N; i++) {
        int[] tmp = new int[2] {1,2};
        for (int j = 0; j < tmp.Length; j++) {
          int n = r.Next(tmp.Length);
          int x = tmp[j];
          tmp[j] = tmp[n];
          tmp[n] = x;
        }

        StringBuilder s = new StringBuilder();
        for (int k = 0; k < tmp.Length; k++) {
          s.Append(tmp[k]);
        }
        
        if(!map.ContainsKey(s.ToString())) {
          map.Add(s.ToString(), 0);
        }
        map[s.ToString()] += 1;
      }

      for(int i = 0; i < map.Count; i++) {
        Console.WriteLine(map.ElementAt(i).Key + " " + map.ElementAt(i).Value);
      }
    }  

    static void fisher(int N) {      
      Dictionary<string, int> map = new Dictionary<string, int>();

      for (int i = 0; i < N; i++) {
        int[] tmp = new int[2] {1,2};
        for (int j = tmp.Length -1; j > 0; j--) {
          int n = r.Next(j+1);
          int x = tmp[j];
          tmp[j] = tmp[n];
          tmp[n] = x;
        }

        StringBuilder s = new StringBuilder();
        for (int k = 0; k < tmp.Length; k++) {
          s.Append(tmp[k]);
        }
        
        if(!map.ContainsKey(s.ToString())) {
          map.Add(s.ToString(), 0);
        }
        map[s.ToString()] += 1;
      }

      for(int i = 0; i < map.Count; i++) {
        Console.WriteLine(map.ElementAt(i).Key + " " + map.ElementAt(i).Value);
      }
    }  
  }
}
