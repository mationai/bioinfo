from superpipe import pipes
import sys
sys.path.append('src')
from helpers import *
from seqlib import *

@pipes
class chap1:
   # Find Indices of min skews of sequence
   def _7():
      dna = 'TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'
      minSkewIndices(dna) >> sepWithSp >> testVs('11 24')

      dna = readlines('dataset_33716_7a.txt')[0]
      minSkewIndices(dna) >> sepWithSp >> testVs('9298 9299')

      dna = readlines('dataset_33716_7.txt')[0]
      minSkewIndices(dna) >> sepWithSp >> testVs('84124')

   # Find hamming distance
   def _8():
      abs(hamming('GGGCCGTTGGT', 'GGACCGTTGAC')) >> testVs(3)

      [s1, s2] = readlines('dataset_33716_8a.txt')
      hamming(s1, s2) >> testVs(751)

      [s1, s2] = readlines('dataset_33716_8.txt')
      hamming(s1, s2) >> testVs(853)

   # Find Indices of ptn in seq with max of d Hamming dist diff 
   def _9():
      dna = 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT'
      indexesIn(dna, 'ATTCTGGA', 3) >> sepWithSp >> testVs('6 7 26 27')

      [ptn, s, n] = readlines('dataset_33716_9a.txt')
      ints = indexesIn(s, ptn, int(n))
      [ints[:4], ints[-4:]]\
         >> testVs([[6, 10, 11, 12], [16389, 16390, 16393, 16394]])

      [ptn, s, n] = readlines('dataset_33716_9.txt')
      ints = indexesIn(s, ptn, int(n))
      # print(ints >> sepWithSp)
      [ints[:4], ints[-4:]]\
         >> testVs([[9, 16, 32, 59], [19525, 19536, 19544, 19559]])

   # Find Hamming dist btw Pattern and every k-mer substring of seq
   def _10():
      indexesIn('TTTAGAGCCTTCAGAGG', 'GAGG', 2) >> len >> testVs(4)

      [ptn, s, n] = readlines('dataset_33716_10a.txt')
      indexesIn(s, ptn, int(n)) >> len >> testVs(12)

      [ptn, s, n] = readlines('dataset_33716_10.txt')
      indexesIn(s, ptn, int(n)) >> len >> testVs(26)

   # Find most frequent k-mers by d dist of mismatch in seq
   def _11():
      [kmrs, cnt] = mostFreqKmers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 1)
      print('cnt:', cnt)
      kmrs >> sepWithSp >> testVs('ATGC ATGT GATG')

      [s, nsStr] = readlines('dataset_33716_11a.txt')
      [k, d] = parseInts(nsStr)
      [kmrs, cnt] = mostFreqKmers(s, k, d)
      print('cnt:', cnt)
      kmrs >> sepWithSp >> testVs('CCCCCC')

      [s, nsStr] = readlines('dataset_33716_11b.txt')
      [k, d] = parseInts(nsStr)
      [kmrs, cnt] = mostFreqKmers(s, k, d)
      print('cnt:', cnt)
      kmrs >> sepWithSp >> testVs('GCACACAGAC GCGCACACAC')

      [s, nsStr] = readlines('dataset_33716_11.txt')
      [k, d] = parseInts(nsStr)
      [kmrs, cnt] = mostFreqKmers(s, k, d)
      print('cnt:', cnt)
      kmrs >> sepWithSp >> testVs('AAAATT TAAAAT')

   # Find most frequent k-mers by d dist of mismatch, account for reverse complement too
   def _12():
      [s, nsStr] = readlines('dataset_33716_12a.txt')
      [k, d] = parseInts(nsStr)
      [kmrs, cnt] = mostFreqKmers(s, k, d, True)
      print('cnt:', cnt)
      kmrs >> sepWithSp >> testVs('CGCGC GCGCG')

      [s, nsStr] = readlines('dataset_33716_12.txt')
      [k, d] = parseInts(nsStr)
      [kmrs, cnt] = mostFreqKmers(s, k, d, True)
      print('cnt:', cnt)
      kmrs >> sepWithSp >> testVs('ATATATA ATCTAGA TATATAT TCTAGAT')

   # Find d-neighbors -- permutations of ptn differ by d or less
   def _13():
      permulate('ACG', 1, True) >> sepWithSp\
         >> testVs('AAG ACA ACC ACG ACT AGG ATG CCG GCG TCG')

      permulate('ACGT', 2, True) >> sepWithSp\
         >> testVs('AAAT AACT AAGA AAGC AAGG AAGT AATT ACAA ACAC ACAG ACAT ACCA ACCC ACCG ACCT ACGA ACGC ACGG ACGT ACTA ACTC ACTG ACTT AGAT AGCT AGGA AGGC AGGG AGGT AGTT ATAT ATCT ATGA ATGC ATGG ATGT ATTT CAGT CCAT CCCT CCGA CCGC CCGG CCGT CCTT CGGT CTGT GAGT GCAT GCCT GCGA GCGC GCGG GCGT GCTT GGGT GTGT TAGT TCAT TCCT TCGA TCGC TCGG TCGT TCTT TGGT TTGT')

      [ptn, n] = readlines('dataset_3014_4.txt')
      S = permulate(ptn, int(n), True)
      [S[:3], S[-3:]] >> testVs([
         parseSeqsStr('AAACGCTGAA AACAGCTGAA AACCACTGAA'),
         parseSeqsStr('TTCTGCTGAA TTGCGCTGAA TTTCGCTGAA')])

      [ptn, n] = readlines('dataset_33716_13.txt')
      S = permulate(ptn, int(n), True)
      # print(S >> sepWithSp)
      [S[:3], S[-3:]] >> testVs([
         parseSeqsStr('AAAAATCATCGT AAAACTCATCGT AAAAGACATCGT'),
         parseSeqsStr('TTCTTTCATCGT TTGTGTCATCGT TTTTGTCATCGT')])

   # Find skew (running sum of ACTG->+-1/0 mapping) of sequence
   def extra():
      getSkew('CATGGGCATCGGCCATACGCC') >> sepWithSp\
         >> testVs('0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2')
      getSkew('GAGCCACCGCGATA') >> sepWithSp\
         >> testVs('0 1 1 2 1 0 0 -1 -2 -1 -2 -1 -1 -1 -1')

   # Find Indices of min skews of sequence in file
   def extra2():
      minSkewInFile('Salmonella_enterica.txt') >> sepWithSp\
         >> testVs('3764856 3764858')

if run(7):
   chap1._7()
if run(8):
   chap1._8()
if run(9):
   chap1._9()
if run(10):
   chap1._10()
if run(11):
   chap1._11()
if run(12):
   chap1._12()
if run(13):
   chap1._13()
   chap1.extra()
   t0 = timeit.default_timer()
   chap1.extra2()
   printRuntime(t0)
   # mbp: ~2s