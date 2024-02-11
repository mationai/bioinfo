from superpipe import pipes
import sys
sys.path.append('src')
from dataexpected import expected
from helpers import *
from seqlib import *
from motifslib import *
from extensions import *

@pipes
class chap2:
   def _1():
      # Expected # of occurrences of a 9-mer in 500 random 1000 long DNA strings is?
      k9Likelihood = 1 / (4 ** 9)
      chances = 500 * (1000 - 9 + 1)
      k9Likelihood * chances >> testVs(1.89208984375)

      S = parseSeqsStr('TTACCTTAAC GATATCTGTC ACGGCGTTCG CCCTAAAGAG CGTCAGAGGT')
      distSum(S, 'AAA') >> testVs(5)
      [ptn, s] = readlines('dataset_5164_1.txt')
      S = parseSeqsStr(s)
      distSum(S, ptn) >> testVs(68)

   # Brute force Motif Enumeration
   def _2():
      S = parseSeqsStr('ATTTGGC TGCCTTA CGGTATC GAAAATT')
      motifsEnum(S, 3, 1) >> testVs(parseSeqsStr('ATA ATT GTT TTT'))

      l = readlines('dataset_33718_2a.txt')
      [k, d] = parseInts(l[0])
      motifsEnum(l[1:], k, d) >> sepWithSp\
         >> testVs('TCAAA TCACA TCAGA TCATA')

      l = readlines('dataset_33718_2b.txt')
      [k, d] = parseInts(l[0])
      S = parseSeqsStr(l[1])
      motifsEnum(S, k, d) >> sepWithSp\
         >> testVs('ACGGC CCCCG CCTCG CGAAC CGCAC CGGAC CGTAC CGTCC CTCGT GCGTA GTAGG TCGGA')

      l = readlines('dataset_33718_2.txt')
      [k, d] = parseInts(l[0])
      S = parseSeqsStr(l[1])
      motifsEnum(S, k, d) >> sepWithSp\
         >> testVs('AACCG ACCCG AGCCG ATCCA ATCCC ATCCG ATCCT ATCTG TCCGG TTCCG')

   def _3exercises():
      maxPossibleMotifScore(10, 15) >> testVs(105)
      maxPossibleMotifScore(10, 12) >> testVs(84)

      S = parseSeqsStr('TCGGGGGTTTTT CCGGTGACTTAC ACGGGGATTTTC TTGGGGACTTTT AAGGGGACTTCC TTGGGGACTTCC TCGGGGATTCAT TCGGGGATTCCT TAGGGGAACTAC TCGGGTATAACC')
      motifsEntropy(S) >> testVs(17.11781961713541) # app marks it wrong

   def _3():
      S = parseSeqsStr('AAATTGACGCAT GACGACCACGTT CGTCAGCGCCTG GCTGAGCACCGG AGTTCGGGACAG')
      [kmr, d] = medianKmer(S, 3)
      print('dbg1 d', d) # 2, 1 for each of the 3rd,4th seqs
      kmr >> testVs('GAC') # debug#1 & section 2.4's last ptn is AGTACGGGACAG, so answer -> ACG

      S = parseSeqsStr('ATA ACA AGA AAT AAC')
      [kmr, d] = medianKmer(S, 3)
      print('dbg3 d', d) # 5, 1 for each seq
      kmr >> testVs('AAA')

      l = readlines('dataset_33718_3b.txt')
      S = parseSeqsStr(l[1])
      [kmr, d] = medianKmer(S, int(l[0]))
      print('dbg5 d', d) # 6
      kmr >> testVs('CGGCGA')

      l = readlines('dataset_33718_3.txt')
      S = parseSeqStrs(l[1:])
      [kmr, d] = medianKmer(S, int(l[0]))
      print('d', d)
      kmr >> testVs('TAAACC') # "wrong" but think it's right

   def _4exercises():
      ptnProbability('TCGTGGATTTCC', {
         'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
         'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
         'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
         'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4],
      }) >> testVs(0.0)
      # 0 since it's product([...]), so if any is 0, product is 0

   # Find most probable k-mer in a sequence
   def _4():
      mostProbableKmer('ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT', 5, {
         'A': [0.2, 0.2, 0.3, 0.2, 0.3],
         'C': [0.4, 0.3, 0.1, 0.5, 0.1],
         'G': [0.3, 0.3, 0.5, 0.2, 0.4],
         'T': [0.1, 0.2, 0.1, 0.1, 0.2]
      }) >> testVs('CCGAG')

      l = readlines('dataset_33718_4a.txt')
      profile = parseProfile(l[2:])
      print('profile', profile)
      mostProbableKmer(l[0], int(l[1]), profile) >> testVs('CGCGGCAGACCAGGA')

      l = readlines('dataset_33718_4.txt')
      profile = parseProfile(l[2:])
      print('profile', profile)
      mostProbableKmer(l[0], int(l[1]), profile) >> testVs('GATTATCTTGAGAC')

   # Greey motif search, no pseudocounts
   def _5():
      dna = parseSeqsStr('GGCGTTCAGGCA AAGAATCAGTCA CAAGGAGTTCGC CACGTCAATCAC CAATAATATTCG')
      greedyMotifSearch(dna, 3, 0) >> testVs(parseSeqsStr('CAG CAG CAA CAA CAA'))

      l = readlines('dataset_33718_5.txt')
      [k, _] = parseInts(l[0])
      seqs = parseSeqsStr(l[1])
      res = greedyMotifSearch(seqs, k, 0) >> sepWithNL # expect NLs
      print('5. res', res)
      res >> testVs(expected._2_5)    

   # Greey motif search
   def _6():
      dna = parseSeqsStr('GGCGTTCAGGCA AAGAATCAGTCA CAAGGAGTTCGC CACGTCAATCAC CAATAATATTCG')
      greedyMotifSearch(dna, 3, 1) >> testVs(parseSeqsStr('TTC ATC TTC ATC TTC'))

      l = readlines('dataset_33718_6.txt')
      [k, _] = parseInts(l[0])
      dna = parseSeqsStr(l[1])
      res = greedyMotifSearch(dna, k, 1) >> sepWithNL # expects NL
      print('6. res', res)
      res >> testVs(expected._2_6)

   def _7():
      seqs = parseSeqsStr('CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACCAA TCCACCAGCTCCACGTGCAATGTTGGCCTA')
      exp = parseSeqsStr('TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG')
      iterRandMotifSearch(seqs, 8) >> testColl(exp)
      gibbsSampler(seqs, 8) >> testColl(exp)

      l = readlines('dataset_33718_7.txt')
      [k, _] = parseInts(l[0])
      exp = parseSeqsStr('ACCGACAGATTGGGT TCCATGCGATTAGGT ACCCATCGATTAGAG TCCCATCGAGCCGGT TCCCACACATTAGGT TCCGTCCGATTAGGT TCCCATCGCCCAGGT TCCCGGAGATTAGGT TCCCCGGGATTAGGT TCCCATCGATTTCTT TCCCATCAGCTAGGT AAACATCGATTAGGT TCCCATCGATGTCGT TCTATTCGATTAGGT TCCCATGAGTTAGGT TCCCAGATATTAGGT TGTGATCGATTAGGT CACCATCGATTAGGC TCCCATATCTTAGGT TCCCATCGATTAACG')
      res = gibbsSampler(parseSeqsStr(l[1]), k)
      print('7. iter res', res >> sepWithSp)
      iterRandMotifSearch(res, k) >> testColl(exp)
      # gibbs same res

   def _7exercises():
      k = 15
      n = 10 # num seqs 
      l = 600 # seq len

      # Probability of capturing a 15-mer in n l-length nucleotides
      p1 = (l - k) / (l+1 - k) # probability of not capturing the kmer in one string
      pr = p1 ** n # p1 * p2 ... * p10, where all probabilities are same
      res = 1.0 - pr
      res >> testVs(0.01693439692911025)

      # Probability of capturing 2 15-mer in n l-length nucleotides
      res = 1.0 - p1**n - (1.0/(l+1-k) * p1**(n-1) * n)
      res >> testVs(0.0001298567056762373)

   def _8():
      # seqs = parseSeqsStr('CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC AATCCACCAGCTCCACGTGCAATGTTGGCCTA')
      # exp = parseSeqsStr('AACGGCCA AAGTGCCA TAGTACCG AAGTTTCA ACGTGCAA')
      # gibbsSampler(seqs, 8, 100) >> testColl(exp) # 100 suggested nor 1k iters work

      l = readlines('dataset_33718_8.txt')
      [k, _, iters] = parseInts(l[0])
      exp = parseSeqsStr('GTAAATAGTCTAGAT GTTCCTCTGCTTGAT GTTAAGGACCTTGAT GTCTCGCTGCTTGAT GTTAAGCTGCTTTGG GTTAAGCGCTTTGAT GTTAAAGAGCTTGAT GTTAGCATGCTTGAT AAGAAGCTGCTTGAT GTTAAGCTGGGCGAT GTTTTACTGCTTGAT GAGGAGCTGCTTGAT GTTATCTTGCTTGAT GTTAAGCTCAATGAT GTTAAGAAACTTGAT CTTAAGCTGCTTGGG GTTAAGCTGCTAACT GTTAATGAGCTTGAT ACTAAGCTGCTTGAA GTTAAGCTGCAAAAT')
      res = gibbsSampler(parseSeqsStr(l[1]), k)
      print('8. gibs res', res >> sepWithSp)
      gibbsSampler(res, k, iters) >> testColl(exp)

if run(1):
   chap2._1()
if run(2):
   chap2._2()
if run(3):
   # chap2._3exercises()
   chap2._3()
if run(4):
   chap2._4exercises()
   chap2._4()
if run(5):
   chap2._5()
if run(6):
   chap2._6()
if run(7):
   chap2._7()
   chap2._7exercises()
if run(8):
   chap2._8()