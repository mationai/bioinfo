from extensions import *
from superpipe import pipes
from pipe import permutations
from functional import seq
from helpers import dPath

A = 'A'
C = 'C'
G = 'G'
T = 'T'
ACGT = [A, C, G, T]
ACGTStr = 'ACGT'

# NOTEs:
# list(gen()) gives "bus error": compiler bug https://github.com/seq-lang/seq/issues/102


def insert[T](a:list[T], i:int, x:T):
   # insert Bug: https://github.com/seq-lang/seq/issues/106 
   if i >= len(a):
      a.append(x)
   else:
      a.insert(i, x)

def slide(dna:str, ptn:str|int, step=1):
   k = len(ptn) if type(ptn) is str else ptn
   return seq(list(dna)).sliding(k, step)\
      .map(lambda l: ''.join(l))\
      .filter(lambda l: len(l) == k)

def hamming(s1:str, s2:str) -> int:
   """ Returns hamming distance between s1 and s2
   """
   return sum(1 for i in range(len(s1)) if s1[i] != s2[i])

@pipes
def indexesSeq(dna:str, ptn:str, maxD=0):
   """ Generates indexes of ptn in dna, with optional max hamming distance d
   """
   for i, s in slide(dna, ptn) >> enumerate:
      match = s == ptn if maxD == 0 else hamming(s, ptn) <= maxD
      if match:
         yield i

## fixme what's d, passed to indexesSeq but not used there
@pipes
def indexesIn(dna:str, ptn:str, maxD=0) -> list[int]:
   """ Like indexesSeq(), but returns actual list, not a generator
   """
   return indexesSeq(dna, ptn, maxD) >> list()

def matched(dna:str, ptn:str, d=0) -> bool:
   """ Returns if ptn is in dna, with optional max hamming distance d
   """
   return any((hamming(s, ptn) <= d for s in slide(dna, ptn)))

def matchedAll(seqs:Seqs, ptn:str, d=0) -> bool:
   """ Returns if (bool) ptn is in every seq in seqs, with optional max hamming distance d
   """
   return all(matched(s, ptn, d) for s in seqs)

def frequentKmers(dna:str, k:int):
   """ Returns most freq k-mers in dna
   """
   cnt = dict[str, int]()
   for kmr in slide(dna, k):
      cnt[kmr] = cnt[kmr]+1 if kmr in cnt else 0

   maxCnt = max(cnt.values())
   return [kmr for kmr in cnt.keys() if cnt[kmr] == maxCnt]

def notCs(c:str) -> Seqs:
   if c == A:
      return [C, G, T]
   elif c == C:
      return [A, G, T]
   elif c == G:
      return [A, C, T]
   elif c == T:
      return [A, C, G]
   return []

def nots(ptn:str) -> Seqs:
   seqs = [notCs(c) for c in ptn]
   return [str(''.join(s)) for s in seqs]

def permulate(ptn:str, d:int, sort=False) -> Seqs:
   # print('ptn',ptn)
   """ 
   d-distance neighbors - all ptns where they differs ptn by d or less
   Eg. 'AC', 1:
     [CC, GC, TC, AA, AG, AT]
   Eg. 'ACT', 1:
     ACC ACT CCT TCT ACG AAT GCT ATT AGT ACA
   Eg. 'ACGT', 1:
     ATGT ACTT ACGG GCGT CCGT AAGT ACGT ACCT AGGT ACAT ACGA TCGT ACGC
   """ 
   if d == 0:
      return [ptn]       
   elif len(ptn) == 1:
      return ACGT

   res = set[str]()
   if d > 1:
      for i, toNot in enumerate(ptn):
         end = '' if i + 1 == len(ptn) else ptn[i +1:]
         _ptn = permulate(ptn[:i] + end, d -1)
         for s in notCs(toNot):
            for p in _ptn:
               # print('add', p[:i] , s , p[i:])
               res.add(p[:i] + s + p[i:])
   else:
      for i, toNot in enumerate(ptn):
         end = '' if i + 1 == len(ptn) else ptn[i +1:]
         for s in notCs(toNot):
            res.add(ptn[:i] + s + end)

   for p in permulate(ptn, d -1):
      res.add(p)
   out = sorted(list(res)) if sort else list(res)
   return out

def reverseComp(dna:str):
   """ Returns reverse complement of dna
   """
   return dna[::-1]\
      .replace('A', 't').replace('T', 'a')\
      .replace('C', 'g').replace('G', 'c').upper()

def mostFreqKmers(dna:str, k:int, d=0, revComp=False,
                  sort=True) -> list[Seqs, int]:
   """ Returns the most most freq k-mers in dna by d-distance or less.
   If revComp=True, account for reverse complement too.
   """
   # cnts = {}
   cnts = dict[str, int]()
   kmrs = slide(dna, k)
   if d <= 0:
      for kmr in kmrs:
         cnts[kmr] = cnts[kmr] +1 if kmr in cnts else 0
   else:
      for p1 in kmrs:
         for p2 in permulate(p1, d):
            cnts[p2] = cnts[p2] +1 if p2 in cnts else 1
         if revComp:
          for p2 in permulate(reverseComp(p1), d):
            cnts[p2] = cnts[p2] +1 if p2 in cnts else 1
   maxCnt = 0 if not cnts else max(cnts.values())
   # print('cnts',cnts)
   res = [kmr for kmr in cnts.keys() if cnts[kmr] == maxCnt]
   out = sorted(res) if sort else res
   return [out, maxCnt]

def permsOfLen(k:int, ptnStr=ACGTStr) -> Seqs:
   """ Permutations of ptnStr of length k
   """
   res = []
   for r in range(len(ptnStr) ** k):
      s = ''
      for j in range(k):
         mod = len(ptnStr) ** (k-j)
         div = len(ptnStr) ** (k-j-1)
         s += ptnStr[(r % mod) // div]
      res.append(str(s))
   return res if k > 0 else Seqs()

def kmerEnums(seqs:Seqs, k:int, d=1, sort=True) -> Seqs:
   """
   All k-mers that appear in seqs with at most d mismatches. 
   Eg. k=3, d=1, seqs=[ATTTGGC TGCCTTA CGGTATC GAAAATT]
   gives [ATA ATT GTT TTT]
   """
   # print('seqs',seqs)
   kmers = set[str]()
   for i, s in enumerate(seqs):
      # print('to slide',s)
      for kmer in slide(s, k):
         for kmr in permulate(kmer, d):
            if len(kmr) == k and\
             matchedAll(seqs[:i] + seqs[i+1:], kmr, d):
            # if matchedAll(seqs[:i] + seqs[i+1:], kmr, d):
               kmers.add(kmr)
   res = list(kmers)
   return sorted(res) if sort else res

def distIn(seqs:Seqs, pattern:str) -> int:
   """ Eg. pattern of AAA:
   d(AAA, TTACCT[TAA]C) = 1
   d(AAA, [ACG]GCGTTCG) = 2
   d(AAA, CCCT[AAA]GAG) = 0 = 1+2 = 3
   """
   k = len(pattern)
   dist = 0
   for s in seqs:
      hamD = len(s)
      for ptn in slide(s, k):
         _hamD = hamming(pattern, ptn)
         if hamD > _hamD:
            hamD = _hamD
      dist += hamD
   return dist

def medianKmer(seqs:Seqs, k:int) -> str:
   """
   Find kmer that w/ least dist over all k-mers in seqs
   Eg.    AAA dist 
   ttaccttAAC 1
   gATAtctgtc 1
   ACGgcgttcg 2
   ccctAAAgag 0
   cgtcAGAggt 1
   dist sum = 5
   """
   dist = len(seqs) * k
   res = ''
   for kmer in permsOfLen(k):
      dist_ = distIn(seqs, kmer)
      if dist > dist_:
         dist = dist_
         res = kmer
   return res

@pipes
def getClumps(dna:str, k:int, L:int, minTimes:int, sort=True) -> Seqs:
   """ All distinct k-mers forming (L, t)-clumps in Genome.
   i.e. distinct k-mers that occurs >= minTimes within L-length in dna
   """ 
   idxs = dict[str, list[int]]()
   for i, kmr in slide(dna, k) >> enumerate:
      if kmr in idxs and len(idxs[kmr]) < minTimes:
         for idx in idxs[kmr]:
            if i - idx + k > L:
               idxs[kmr].pop(0)
         idxs[kmr].append(i)
      elif kmr not in idxs:
         idxs[kmr] = [i]
   res = [kmr for kmr in idxs.keys() if len(idxs[kmr]) >= minTimes]
   return sorted(res) if sort else res

@pipes
def getSkew(dna:str, skew0=0) -> list[int]:
   """ Returns skew (running sum of ACTG->+-1/0 mapping) of sequence
   """
   def nextSkew(k, n=0):
      d = {
         'A': 0,
         'T': 0,
         'C': -1,
         'G': 1,
      }
      return d[k] + n
   skew = [skew0]
   for i, k in list(dna) >> enumerate:
      skew.append(nextSkew(k, skew[i]))
   return skew

def minSkewIndices(dna:str):
   """ Returns indices of min skews of sequence
   """
   low = min(getSkew(dna))
   return [i for i, n in enumerate(getSkew(dna)) if n == low]

def minSkewInFile(filename:str) -> list[int]:
   """ Returns indices of min skews of sequence in file
   FASTA('path/Salmonella_enterica.fa') # will work in later releases.
   """
   f = open(dPath(filename), 'r')
   skews = list[list[int]]()
   lows = list[int]()
   low = 0
   lineLen = 0 

   for i, line in enumerate(f):
      if i == 0:
         continue
      s = line.strip('\r\n')
      if i == 1:
         lineLen = len(s)
      skew = getSkew(s, 0 if i==1 else skews[-1][-1])
      low = min(low, (min(skew)))
      skews.append(skew)

   for i, skew in enumerate(skews):
      for j, n in enumerate(skew):
         if n == low:
            lows.append(i * lineLen + j)
   return lows