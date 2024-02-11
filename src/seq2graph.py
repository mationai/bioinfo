from extensions import *
from seqlib import slide
from helpers import adjGraphStr


class Seq2Graph:
   sg: Graph
   def __init__(self, dna:Seqs, method:str, k:int):
      """ Returns de Bruijn graph via:
      """
      match method:
         case 'seq':
            self.sg = self.fromSeq(dna[0], k)
         case 'overlap': # Don't use bad result, use 'kmers'
            self.sg = self.viaOverlap(dna)
         case _:
            self.sg = self.fromKmers(dna)

   def fromKmers(self, kmers:Seqs) -> Graph:
      """ Returns Graph Graph where each node value is of length kmer-length - 1 
      Eg. GAGG CAGG GGGG GGGA CAGG AGGG GGAG =>
      AGG -> GGG
      CAG -> AGG,AGG
      GAG -> AGG
      GGA -> GAG
      GGG -> GGA,GGG
      Can't handle 2-mers, eg. AG AT AA GA GG GT TA TG TT AT -> 
      AGATAATGGTT (overlap worst, -> AATGGAAGGTGATATTGTAGAGTTAA)
      AAGTTGGATAT - expected
      """
      if len(kmers) == 0:
         raise ValueError('Can not build adjacencyGraph from empty seqs')

      # Set k to be length of k-mer to look at, not length input kmers 
      k = len(kmers[0]) - 1
      if k < 1:
         raise ValueError('Can not build adjacencyGraph via "kmers" method with k=0')
      res = Graph()
      for kmer in kmers:
         x, y = [str(s) for s in slide(kmer, k)]
         if x in res:
            res[x].append(y)
         else:
            res[x] = [y]
      return res

   def viaOverlap(self, dna:Seqs) -> Graph:
      """ Eg: ATGCG GCATG CATGC AGGCA GGCAT GGCAC =>
      IndexError: list index out of range
      GCATG -> CATGC
      CATGC -> ATGCG
      GGCAT -> GCATG
      AGGCA -> GGCAT,GGCAC
      """
      if len(dna) < 2:
         raise ValueError('Can not build adjacencyGraph from seqs of length 0 - 1')
      seqs = [str(s) for s in set(dna)]
      res = Graph()
      for i, x in enumerate(seqs):
         for y in seqs[0:i] + seqs[i+1:]:
            if x[1:] == y[:-1]:
               if x not in res:
                  res[x] = Strs()
               res[x].append(y)
      return res

   def fromSeq(self, dna:str, k=4) -> Graph:
      """ Eg: CTTCTTC, 4 =>
      CTT -> TTC,TTC
      TTC -> TCT
      TCT -> CTT
      """
      if k < 2:
         raise ValueError('Can not build adjacencyGraph from k of 0 - 1')
      seqs = [str(s) for s in slide(dna, k-1)] # result are edge labels, so k-1
      res = Graph()
      for i, x in enumerate(seqs[:-1]):
         if x not in res:
            res[x] = Strs()
         res[x].append(seqs[i+1])
      return res

   def toStrs(self) -> [str]:
      return [f"{x}: {' '.join(self.sg[x])}" for x in self.sg]

   def toStr(self) -> str:
      return adjGraphStr(self.sg)

   def __str__(self) -> str:
      return adjGraphStr(self.sg)
