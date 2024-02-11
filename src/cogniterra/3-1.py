from superpipe import pipes
import sys
sys.path.append('src')
from helpers import *
from graphfns import *
from seqlib import *
from extensions import *
from seq2graph import Seq2Graph

@pipes
class chap3:
   # Genome Path Problem
   def _2():
      S = [kmr for kmr in slide('CAATCCAAC', 5)]
      S >> testVs(parseSeqsStr('CAATC AATCC ATCCA TCCAA CCAAC'))

      [k, seq] = readlines('dataset_34403_2.txt')
      S = [kmr for kmr in slide(seq, 100)]
      # print(S >> sepWithSp)
      S[0] >> testVs('GACACAGAATTGGACCTCTTGGATCGGTCGTACCTGATTAAGCGCGGAAATTTCGTGGGTCTAAGCGCGAGGTATGTCGTTACCAAGCCAGCATGGTCAG')
      S[-1] >> testVs('AGGATAGGTTGAACCCGGTGTTAGAGCAGCAATAGGGAACGCTCCAGGGTACGCGGGTCAAGTTGGACTCCACCGAGCCTTCCGTCTATACACCCTCTTC')

   def _3():
      s = parseSeqsStr('ACCGA CCGAA CGAAG GAAGC AAGCT')
      mergeOrderedSeqs(s) >> testVs('ACCGAAGCT')

      s = parseSeqsStr('CAATC AATCC ATCCA TCCAA CCAAC')
      mergeOrderedSeqs(s) >> testVs('CAATCCAAC')

      s0 = readlines('dataset_34403_3.txt')[0]
      s = parseSeqsStr(s0)
      S = mergeOrderedSeqs(s)
      S[:45] >> testVs('GCCAAGTGATCGCATGGGAAGGGCACATTGCGAGAGCGACCCATA')

   def _4():
      s = parseSeqsStr('ATGCG GCATG CATGC AGGCA GGCAT GGCAC')
      r = Seq2Graph(s, 'overlap', 0)
      r.toStrs() >> testColl([ # sometimes fail
         'GCATG: CATGC',
         'CATGC: ATGCG',
         'GGCAT: GCATG',
         'AGGCA: GGCAC GGCAT'])

      s = parseSeqsStr('CT TT TT')
      r = Seq2Graph(s, 'overlap', 0)
      r.toStr() >> testVs('CT: TT')

      s = readlines('dataset_34403_4a.txt')
      r = Seq2Graph(s, 'overlap', 0)
      res = r.toStrs()
      testVs(len(res), 4975)
      testIn(res, 'AGTCTGAACTTAGTGCTGCCTCAGA: GTCTGAACTTAGTGCTGCCTCAGAT')
      testIn(res, 'TGGCAGCAGGCATCTTGGCCAGTAT: GGCAGCAGGCATCTTGGCCAGTATC')

      s0 = readlines('dataset_34403_4.txt')[0]
      s = parseSeqsStr(s0)
      r = Seq2Graph(s, 'overlap', 0)
      res = r.toStrs()
      # print(res)
      testVs(len(res), 4975)
      testIn(res, 'ACTATGGATGTTTCGTCAACAGGAG: CTATGGATGTTTCGTCAACAGGAGC')
      testIn(res, 'CTTTAGTTTCTCGTTAAGAGATTCT: TTTAGTTTCTCGTTAAGAGATTCTA')

   def _5():
      r = Seq2Graph(['AAGATTCTCTAAGA'], 'seq', 4)
      r.toStrs() >> testColl([ # sometimes fail
         'AAG: AGA AGA',
         'AGA: GAT',
         'ATT: TTC',
         'CTA: TAA',
         'CTC: TCT',
         'GAT: ATT',
         'TAA: AAG',
         'TCT: CTC CTA',
         'TTC: TCT',
      ])

      r = Seq2Graph(['GCTTCTTC'], 'seq', 4)
      r.toStrs() >> testColl([
         'GCT: CTT',
         'TTC: TCT',
         'CTT: TTC TTC',
         'TCT: CTT'])

      [k, s] = readlines('dataset_34403_5a.txt')
      r = Seq2Graph([s], 'seq', int(k))
      res = r.toStrs()
      testVs(len(res), 1988)
      testIn(res, 'TGACATAACGT: GACATAACGTT')
      testIn(res, 'GCACGAGCCAG: CACGAGCCAGA')

      [k, s] = readlines('dataset_34403_5.txt')
      r = Seq2Graph([s], 'seq', int(k))
      # print(r)
      res = r.toStrs()
      testVs(len(res), 1988)
      testIn(res, 'TCTCAGACGGA: CTCAGACGGAG')
      testIn(res, 'TACCATATGCT: ACCATATGCTA')

   def _6():
      s = parseSeqsStr('GAGG CAGG GGGG GGGA CAGG AGGG GGAG')
      r = Seq2Graph(s, 'kmer', 0)
      r.toStrs() >> testColl([
         'CAG: AGG AGG',
         'AGG: GGG',
         'GGA: GAG',
         'GAG: AGG',
         'GGG: GGG GGA'])

      dna = readlines('dataset_34403_6a.txt')
      r = Seq2Graph(dna, 'kmer', 0)
      res = r.toStrs()
      testVs(len(res), 2171)
      testIn(res, 'AACTATTGGATATCGAAAC: ACTATTGGATATCGAAACA')
      testIn(res, 'TTTAAACGGCCCGGATGCA: TTAAACGGCCCGGATGCAG')

      s0 = readlines('dataset_34403_6.txt')[0]
      dna = parseSeqsStr(s0)
      r = Seq2Graph(dna, 'kmer', 0)
      # print(r)
      res = r.toStrs()
      testVs(len(res), 2172)
      testIn(res, 'ACATTTGTCCAGACCAACC: CATTTGTCCAGACCAACCA')
      testIn(res, 'AAGATATGTAATATTATAG: AGATATGTAATATTATAGA')

if run(2):
   chap3._2()
if run(3):
   chap3._3()
if run(4):
   chap3._4()
if run(5):
   chap3._5()
if run(6):
   chap3._6()
