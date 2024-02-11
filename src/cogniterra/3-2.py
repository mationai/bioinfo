from superpipe import pipes
import sys
sys.path.append('src')
from extensions import *
from helpers import *
from seq2graph import Seq2Graph
from graphlib import *
from seqlib import permsOfLen

@pipes
class chap3:
   # Find Eulerian Cycle
   def _7():
      s = '''0: 3
      1: 0
      2: 1,6
      3: 2
      4: 2
      5: 4
      6: 5,8
      7: 9
      8: 7
      9: 6'''
      S = eulerianPath(parseGraph(s), 'cycle', '')
      S >> pathArrowStr >> rmSpaces >> testVs('2->1->0->3->2->6->8->7->9->6->5->4->2')

      s = '1: 2|2: 1,2|0: 1'
      S = eulerianPath(parseGraph(s, '|'), 'cycle', '')
      S >> pathArrowStr >> rmSpaces >> testVs('')

      s = '''1: 10
      10: 2,3,4
      2: 1
      3: 10
      4: 5
      5: 10'''
      S = eulerianPath(parseGraph(s), 'cycle', '')
      S >> pathArrowStr >> rmSpaces >> testVs('10->2->1->10->3->10->4->5->10')

      s = '''0: 1,2,3,4
      1: 0,2,3,4
      2: 0,1,3,4
      3: 0,1,2,4
      4: 0,1,2,3'''
      S = eulerianPath(parseGraph(s), 'cycle', '')
      S >> pathArrowStr >> rmSpaces >> testVs('0->1->0->2->0->3->0->4->1->2->1->3->1->4->2->3->2->4->3->4->0')

      s = '0: 3,1|1: 2|2: 0|3: 0'
      S = eulerianPath(parseGraph(s, '|'), 'cycle', '')
      S >> pathArrowStr >> rmSpaces >> testVs('0->3->0->1->2->0')

      gs = readGraphStr('dataset_203_2.txt') 
      o = read('2-2-1.2.2-out.txt')
      eulerianPath(gs, 'cycle', '') >> pathArrowStr >> rmSpaces >> testVs(o)

   # Find Eulerian Path
   def _8():
      s = '''0: 2
      1: 3
      2: 1
      3: 0,4
      6: 3,7
      7: 8
      8: 9
      9: 6'''
      S = eulerianPath(parseGraph(s), 'path', '')
      S >> pathArrowStr >> rmSpaces >> testVs('6->7->8->9->6->3->0->2->1->3->4')

      gs = readGraphStr('dataset_203_6.txt')
      o = read('2-2-1.2.6-out.txt')
      eulerianPath(gs, 'path', '') >> pathArrowStr >> rmSpaces >> testVs(o)

   # Reconstruct Genome from scrambled seqs
   def _9():
      S =  parseSeqsStr('CTTA ACCA TACC GGCT GCTT TTAC')
      genomeFromSeqs(S, 'path', '') >> testVs('GGCTTACCA')

      S = parseSeqsStr('AAC AAC ACG ACT CGA GAA')
      genomeFromSeqs(S, 'path', '') >> testVs('AACGAACT')

      S = parseSeqsStr('CTAC CTCC TCCT ACTC CCTC CCTA TACT')
      genomeFromSeqs(S, 'path', '') >> testIn(['CCTACTCCTC', 'CCTCCTACTC'])

      S = parseSeqsStr('CCC CCC CCC TCC CCC CCG CCC CCC CCC')
      genomeFromSeqs(S, 'path', '') >> testVs('TCCCCCCCCCG')

      S = parseSeqsStr('AG AT AA GA GG GT TA TG TT AT')
      genomeFromSeqs(S, 'path', '') >> testIn(['AAGTTGGATAT', 'AGATAATGGTT'])

      l = readlines('dataset_203_7.txt')
      o = read('2-2-1.2.7-out.txt')
      parseSeqStrs(l[1:]) >> genomeFromSeqs('path', '') >> testVs(o)

   # 1.2.9 - only 1 k-Universal Circular Path and its reverse, so answer = 2

   # Find k-Universal Circular String
   def _10():
      S = permsOfLen(3, '01')
      genomeFromSeqs(S, 'universal', '') >> testVs('00101110')

      S = permsOfLen(4, '01')
      genomeFromSeqs(S, 'universal', '') >> testVs('0010011010111100')

      S = permsOfLen(9, '01')
      genomeFromSeqs(S, 'universal', '') >> testVs('10100000000010000001010000010010000011000000011010000100010000101010000110010000111000000111010100010010100011000100011010100100010100101010100110000100110010100111000100111011000001011000011011000101011000110011000111001000111011100001011100011011100101011100110011100111100000111100100100101100100110100100111101000101101000111101100101101100110101100111110000111110100101110100110110100111111000101111000111111100101111100110111101010101101010111101110101101110110101110111110110110111111010111111111011110011')

   # 1.3 - Assembling Genomes from Read-Pairs

   # Generate (k, d)-mer Composition of a pattern 
   # Reconstruct genome from Pairs
   def _11():
      ptn = 'TAATGCCATGGGATGTT'
      s = pairedComposition(ptn, 3, 1) >> toPairsStr
      s >> testVs('(AAT|CCA) (ATG|CAT) (ATG|GAT) (CAT|GGA) (CCA|GGG) (GCC|TGG) (GGA|GTT) (GGG|TGT) (TAA|GCC) (TGC|ATG) (TGG|ATG)')

      s = pairedComposition(ptn, 3, 2) >> toPairsStr
      s >> testVs('(AAT|CAT) (ATG|ATG) (ATG|ATG) (CAT|GAT) (CCA|GGA) (GCC|GGG) (GGG|GTT) (TAA|CCA) (TGC|TGG) (TGG|TGT)')

      sp = pairsStrToStrPairs('GAGA|TTGA TCGT|GATG CGTG|ATGT TGGT|TGAG GTGA|TGTT GTGG|GTGA TGAG|GTTG GGTC|GAGA GTCG|AGAT', ' ')
      genomeFromPairs(sp, 4, 2) >> testVs('GTGGTCGTGAGATGTTGA')

      sp = pairsStrToStrPairs('TCA|GCA TTC|TGC AAT|CAT ATT|ATG', ' ')
      genomeFromPairs(sp, 3, 1) >> testVs('AATTCATGCA')

      sp = pairsStrToStrPairs('GG|GA GT|AT TG|TA GA|AC AT|CT', ' ')
      genomeFromPairs(sp, 2, 1) >> testVs('GGTGATACT')

      sp = pairsStrToStrPairs('GTTT|ATTT TTTA|TTTG TTAC|TTGT TACG|TGTA ACGT|GTAT CGTT|TATT', ' ')
      genomeFromPairs(sp, 4, 2) >> testVs('TTTACGTTTGTATTT')

      sp = pairsStrToStrPairs('GGG|GGG AGG|GGG GGG|GGT GGG|GGG GGG|GGG', ' ')
      genomeFromPairs(sp, 3, 2) >> testVs('AGGGGGGGGGGT')

   # NOTE: not correct. Using maxNonBranchingPaths should help 
   def _12():
      l = readlines('2-2-1.3.16-testin.txt')
      o = readlines('2-2-1.3.16-testout.txt')
      i = parseInts(l[0])
      sp = pairStrsToStrPairs(l[1:], '\n')
      genomeFromPairs(sp, i[0], i[1]) >> testVerbose(o[0])

      l = readlines('dataset_204_16.txt')
      i = parseInts(l[0])
      sp = pairStrsToStrPairs(l[1:], '\n')
      genomeFromPairs(sp, i[0], i[1]) >> echo

   # Construct genome from ordered pairs or "string spelled by gapped patterns"
   def _13():
      sp = pairsStrToStrPairs('GACC|GCGC ACCG|CGCC CCGA|GCCG CGAG|CCGG GAGC|CGGA', ' ')
      genomeFromOrderedPairs(sp, 4, 2) >> testVs('GACCGAGCGCCGGA')

      l = readlines('dataset_6206_4.txt')
      i = parseInts(l[0])
      sp = pairStrsToStrPairs(l[1:], '\n')
      genomeFromOrderedPairs(sp, i[0], i[1]) >> testVs('AAGCCGTCGCCCAGGAGTTCGATGGCGGCGACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGATGCACACAGGCCGCATGTTTTATAATAATAATTATGGGATTACCTCTTATTACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGAATGAACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGATTGGGTTTGAGTCACGCGGGAGTCGTCCTACGAACACTCTTCATTGTTCGCAGCCTACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGAGCCGCACACTTATCCATTAGCAGAACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGATTCTTGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGAATTAAAACACGTCAGTCGTGGAATTAAGATATTGGTCTCCACGCGACTTATCCATTAGCAACTTATCCATTAGCAGATCTATAGACGCGCGGGACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGAGTCTTAAAGTTTCTTGAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGAGTGGTTGTGTACAATTCGAGTGTTGCTTGACATGGCAGTATCACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGAAGACCAGTTAAGCTGGTCCGTGACGTTCACCATTTCAAGCCATAAATCGACTTCCATGCACCACCAGCGCGCACAGAATCCGCTTCAACTCGTACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGACTCGCAAACCTGACTAATTAAAAGCTGCATCCACAGGATCTGAGACGAATTGTCTCTAAAAGAAGACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGACAATGGACCATAACTAAATTTCACTCGGCCTAAGATGTCACGAATGCTAGACCTCGGGCCTCTCTGGCCACAAATCGTACTCTTGTAACCTATGGCATCATATCCGGAGAGCTCCTGCATTTCAAGAAGGACACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGATTTCTTGACCTTCGACACCTATACCAAATTCCTACCCAATAGATCAAACGTGGCTAAAGTAATGGTGGGTTAACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGAACGTCATTATTATCAGCAGCGACGTCTCTTTGTTCACACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGAATAGTCCGCCGGCACTACACGTCTAGGCGGCAGTTTTCTAAGAATTGTATGGCGACTCAGTCATTGCATCGCGTAGCCTTGACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGACGAACAACCGGTCTTCCTGAATAAGCGTAAGAGTGAGGACACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGAACACGTAAAAGAGTCGGGCTATAACCACTAATTCCGTCTACGAGCTATGGAAACCAAAAAGTTTCTTAGGACTAGCGACCCCCCCTTACCCCGGGTGTGTACAGAAAGGTGCCTTAACGCATCCCAACTTATCCATTAGCAGATCTATAGACGCGCGGGGTCTTAAAGTTTCTTGAAGACCCTACACACCCGGGTCCAGTCGGACATTTCTGGTATTGCCATGACGATCCGCCGACCGGCCAGACCTACGGGCTCTCGGCAGGCTGCACGCACAACACCTCCGCTGCGCGGCCCGAAAGTCAG')

   # Find maximal non-branching paths
   def _14():
      S = parseGraph('1: 2|2: 3|3: 4,5|6: 7|7: 6', '|')
      maxNonBranchingPaths(S) >> pathsArrowStr >> testVs('1: 2: 3\n3: 4\n3: 5\n6: 7: 6')

      S = readGraphStr('dataset_6207_2.txt')
      o = read('2-2-1.6.2-out.txt')
      maxNonBranchingPaths(S) >> pathsArrowStr >> testVs(o)

   # Generate Contig - long, contiguous segments of non-branching genome
   def _15():
      ptn = 'TAATGCCATGGGATGTT' # See 1.4.4
      contigsFromSeq(ptn, 3) >> testColl('TAAT TGTT TGCCAT ATG ATG ATG TGG GGG GGAT'.split(' '))

      kmers = parseSeqsStr('ATG ATG TGT TGG CAT GGA GAT AGA')
      contigs(kmers) >> testColl('AGA ATG ATG CAT GAT TGGA TGT'.split(' '))

      kmers = readseqs('dataset_205_5.txt')
      r = contigs(kmers)
      # r >> sepWithSp >> echo
      r[0] >> testVs('CATCCAAAGGGACGAGATCAAAACCGGTAGTACCAGCCGTTTTCGACACTGGGAATTTCCCTAATAAGA')
      r[-1] >> testVs('TTACAAAACGCTCCCACACAGTATCACCTTGTTGTGCAAAATCGAAGGTCTTGGTCTGGCGCGGCTGTACCAGATTAGCAAGAACGCGCCCAGGATCTTACCCAGATCCCCAGAGGTCGTTCCGTTTGATGCATC')

if run(8):
   chap3._8()
if run(9):
   chap3._9()
if run(10):
   chap3._10()
if run(11):
   chap3._11()
if run(12):
   chap3._12()
if run(13):
   chap3._13()
if run(14):
   chap3._14()
if run(15):
   chap3._15()