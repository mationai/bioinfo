from math import ceil, log2
import random
from extensions import Profile, SeqFloatDict
from helpers import product
from seqlib import *


def acgtHistogram(motifs:Seqs, initCnt=1) -> dict[str, list[int]]:
   """ Returns {
     'A': [position-0 count, position-1 count, ...]
     ...
   }
   """
   motifLen = len(motifs[0])
   res = {
      A: [initCnt for j in range(motifLen)],
      C: [initCnt for j in range(motifLen)],
      G: [initCnt for j in range(motifLen)],
      T: [initCnt for j in range(motifLen)],
   }
   for i in range(len(motifs)):
      for j in range(motifLen):
         c = motifs[i][j]
         res[c][j] += 1
   return res

def motifsEnum(seqs:Seqs, k:int, d=1, sort=True) -> Seqs:
   """
   Find motifs in seqs with at most d mismatch, brute force.
   Eg. k=3, d=1, seqs=[ATTTGGC TGCCTTA CGGTATC GAAAATT]
   gives [ATA ATT GTT TTT]
   """
   # print('seqs',seqs)
   kmers = set[str]()
   for i, seq in enumerate(seqs):
      # print('to slide',s)
      for kmer in slide(seq, k):
         for kmr in permulate(kmer, d):
            if len(kmr) == k and\
             matchedAll(seqs[:i] + seqs[i+1:], kmr, d):
            # if matchedAll(seqs[:i] + seqs[i+1:], kmr, d):
               kmers.add(kmr)
   res = list(kmers)
   return sorted(res) if sort else res

def motifsProfile(motifs:Seqs, initCnt=1) -> Profile:
   """ 
   Motifs profile is a matrix of probabilities of each nucleotide
     at each position in motifs, where probability is:
      count of nucleotide / count of motifs
   Returns {
     'A': [probabilities...]
     ...
   }
   """
   counts = acgtHistogram(motifs, initCnt)
   return {
      c: [x/float(len(motifs)) for x in counts[c]]
         for c in counts
   }

def motifsScore(motifs:Seqs, initCnt=1) -> int:
   """ Returns
   sum of score in each col, where score is count of occurences of non-most common nucleotide in col.
   """
   if len(motifs) == 0:
      raise ValueError('Unexpected motifs length of 0')
   motifLen = len(motifs[0])
   maxcount = [0 for i in range(motifLen)]
   counts = acgtHistogram(motifs, initCnt)

   for i in range(motifLen):
      for c in ACGT:
         if maxcount[i] < counts[c][i]:
            maxcount[i] = counts[c][i] 
   return sum(len(motifs) - maxcount[i] for i in range(motifLen))

def motifsEntropy(motifs:Seqs) -> float:
   """ Returns Entropy = measure of the
   uncertainty of a probability distribution of motifs profile
   """
   res = 0.0
   profile = motifsProfile(motifs)
   for i in range(len(motifs[0])):
      entropy = 0.0
      for c in ACGT:
         p = profile[c][i]
         if p > 0:
            entropy += abs(p * log2(p))
      res += entropy
   return res

def ptnProbability(ptn:str, profile:Profile) -> float:
   """ Probability of ptn given profile
   """
   try:
      return product([profile[c][i] for i, c in enumerate(ptn)])
   except IndexError:
      raise IndexError("profile[A|C|G|T]'s value length < ptn's length")

def kmerProbabilities(dna:str, profile:Profile, k:int) -> list[float]:
   """ All probabilities of kmers in dna
   """
   return [ptnProbability(kmer, profile) for kmer in slide(dna, k)]

def mostProbableKmer(ptn:str, k:int, profile:Profile) -> str:
   """ Most probabile k-mer of given ptn & profile, aka Profile-most probable k-mer
   """
   maxPrb = -1.0
   res = ''
   for kmer in slide(ptn, k):
      prob = ptnProbability(kmer, profile)
      if prob > maxPrb:
         maxPrb = prob
         res = kmer
   return res

def mostProbableSeq(seqs:Seqs, profile:Profile) -> str:
   """ Most probabile k-mer of given ptn & profile, aka Profile-most probable k-mer
   """
   maxPrb = -1.0
   res = ''
   for dna in seqs:
      prob = ptnProbability(dna, profile)
      if prob > maxPrb:
         maxPrb = prob
         res = dna 
   return res

def greedyMotifSearch(seqs:Seqs, k:int, initCnt=1):
   """
   Like medianKmer, but much faster, less accurate.
   Accuracy somewhat mitigated by having non-0s as profile, 
    done via initCnt=1 (pseudocounts).

   Returns a kmer for each str where kmer has best match across seqs. Eg: 
   ggcgttCAGgca 1 
   aagaatCAGtca 1
   CAAggagttcgc 0 
   cacgtCAAtcac 0
   CAAtaatattcg 0
        Score = 2
   http://www.mrgraeme.co.uk/greedy-motif-search/
   https://cogniterra.org/lesson/29870/step/1?unit=21968
   Note: The performance of function deteriorate if first string in Dna does not
    contain k-mers similar to the consensus motif. Thus, biologists often run
    it multiple times, shuffling the strings in Dna each time.
   """
   seqsLen = len(seqs)
   bestMotifs = [seqs[i][:k] for i in range(seqsLen)]

   for kmer in slide(seqs[0], k):
      motifs = [kmer]
      for i in range(1, seqsLen):
         profile = motifsProfile(motifs, initCnt)
         motifs.append( mostProbableKmer(seqs[i], k, profile) )
      if motifsScore(motifs) < motifsScore(bestMotifs, initCnt):
         bestMotifs = motifs 
   return bestMotifs

def bestProfileMotifs(profile:Profile, seqs:Seqs) -> Seqs:
   """ Find the best motif according to a profile and a list of dna strings
   """
   k = len(profile[A])
   res = seqs[:]
   for j in range(len(seqs)):
      score = -1.0
      for kmer in slide(seqs[j], k):
         score_ = ptnProbability(kmer, profile) # probability of occurence
         if score_ > score:
            score = score_
            res[j] = kmer
   return res

def iRandj(seqs:Seqs, k:int):
   seqLen = min([len(s) for s in seqs])
   for i in range(len(seqs)):
      yield i, random.randint(0, seqLen-k)

def _randMotifSearch(seqs:Seqs, k:int) -> Seqs: 
   """ Returns motifs of k-len
   """
   seqLen = len(seqs[0])
   res = motifs = [seqs[i][j:j+k] for i, j in iRandj(seqs, k)]
   while True:
      profile = motifsProfile(motifs)
      motifs = bestProfileMotifs(profile, seqs)
      if motifsScore(motifs) < motifsScore(res):
         res = motifs
      else:
         return res

def iterRandMotifSearch(seqs:Seqs, k:int, iters=1000) -> Seqs:
   """ 
   Note: Often wrong res on small inputs, tried iters=2000, no help

   Returns motifs of k-len, by running randMotifSearch n times
   """
   score = len(seqs[0]) * len(seqs)
   res = []
   for _ in range(iters):
      motifs = _randMotifSearch(seqs, k)
      score_ = motifsScore(motifs)
      if score_ < score:
         score = score_
         res = motifs
   return res

def randomKey(seqProb:SeqFloatDict) -> str:
   """
   Randomly select a key from seqProb Eg.: {'ACGTG': .5, ...}, weighted by its probs.
   """
   aggr = sum(seqProb.values())
   profile = {dna: seqProb[dna]/aggr for dna in seqProb}
   y = 0.0
   profileSums = dict[str, float]()

   for dna in profile: #need to add them together so when a number is randomly selected, there will be defined ranges
      y += profile[dna]
      profileSums[dna] = y
   randf = random.uniform(0, 1)
   for dna in profileSums:
      if randf <= profileSums[dna]:
         return dna

def genKmer(dna:str, profile:Profile, k:int) -> str:
   """
   Generate a k-mer weighted on the probability profile of it occuring in the dna string
   """
   seqProb = {dna[i:i+k]: ptnProbability(dna[i:i+k], profile) for i in range(len(dna)-k+1)}
   return randomKey(seqProb)

def gibbsSampler(seqs:Seqs, k:int, iters=3000, seedIters=300) -> Seqs:
   """
   Note: Often wrong res on small inputs, tried iters=5000, seedIters=500, no help

   GibbsSampler explores just a small subset of solutions, it may "get stuck" in a local optimum.
   For this reason, similarly to randMotifSearch, it should be run many times with the hope that
    one of these runs will produce the best-scoring motifs.
   When we run GibbsSampler 2,000 times on the Subtle Motif Problem with implanted 15-mer
    AAAAAAAAGGGGGGG (each time with new randomly selected k-mers for N = 200 iterations),
    it returns a collection Motifs with consensus string AAAAAAgAGGGGGGt and Score(Motifs)
    equal to 38. This score is even lower than the score of 40 expected from the implanted motifs!
   """
   bestScore = len(seqs[0]) * len(seqs)
   res = []

   # run randMotifSearch few times to seed initial motifs
   for i in range(seedIters):
      motifs = _randMotifSearch(seqs, k)
      score = motifsScore(motifs)
      if score < bestScore:
         bestScore = score
         res = motifs 

   for run in range(iters):
      i = random.randint(0, len(seqs)-1)
      motifs = res[:]
      motifs.pop(i)
      profile = motifsProfile(motifs)
      kmer = genKmer(seqs[i], profile, k)
      insert(motifs, i, kmer)

      if motifsScore(motifs) < motifsScore(res):
         res = motifs
   return res

def maxPossibleMotifScore(n:int, motifLen:int) -> int:
   """
   Maximum score occurs when the motif matrix is "least conserved",
    i.e. when the max # of one letter in a column is minimized.
   This occurs when each letter occurs numMofits/4 times in each col. 
   If n is not a multiple of 4, then there will be one letter that 
    occurs more than n/4 times in a column. 
    This letter occurs ceiling(n/4) times.
   Multiply this by motifLen columns for maximum possible score.
   """
   return int(motifLen * (n - ceil(n/4)))

## aliases ##
def motifEnum(seqs:Seqs, k:int, d:int) -> Seqs:
   return motifsEnum(seqs, k, d)
