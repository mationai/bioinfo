from superpipe import pipes
import sys
sys.path.append('src')
from dataexpected import expected
from extensions import *
from helpers import *
from graphfns import *
from alignment import *

@pipes
class chap5:
   #* Coin Change
   def _2():
      ints = [50,25,20,10,5,1]
      minPickDP(40, ints) >> testVs(2)
      ints = [15,14,5,3,1]
      minPickDP(18347, ints) >> testVs(1224)
      ints = [5,3,1]
      minPickDP(19697, ints) >> testVs(3941)
   
   #* Manhattan Tourist
   def _3():
      mats = parseMats('''1 0 2 4 3
4 6 5 2 1
4 4 5 2 1
5 6 8 5 3
-
3 2 4 0
3 2 4 2
0 7 3 3
3 3 0 2
1 3 2 2''')
      manhattanTourist(mats[0], mats[1]) >> testVs(34)

      mats = parseMats('''4 1 0 4 2
4 3 2 0 2
1 0 3 1 0
1 4 4 2 4
1 0 4 1 3
-
4 4 2 4
4 1 1 4
0 1 4 4
1 3 1 4
4 3 3 3
2 0 2 3''')
      manhattanTourist(mats[0], mats[1]) >> testVs(28)

      l = readlines('dataset_34405_3.txt')
      mats = parseMats('\n'.join(l[1:]))
      manhattanTourist(mats[0], mats[1]) >> testVs(60)

   #* Longest Common Subsequence
   def _4():
      longestCommonSubseq('AACCTTGG', 'ACACTGTGA') >> testVs('AACTTG')

      l = readlines('dataset_34405_4a.txt')
      longestCommonSubseq(l[0], l[1]) >> testVs('GGTACACATCAAGTCCAGGCCAAGTTATTTTCGAATGTGAAGCATGTTTCATACAGACCCACCTCCAGTTACATCGACTTTGGCCTACTTCTGACTGAAATGACTGATAGGAAATAAACTACGTTCCGACATATCGACTTTTTCGAATTTATAGTAGTGGCCATCCTATTTTGTGTTCTGACAATGCGTACATTACATATGACCTTACTTCGAGACGGTGATATTGTATCGAAGGCAAAACCTCCGGTGGCAGCGCGCTGCAATCACGAATCGACCCAGGTCGGTGGAACTGCGCCCTACACTTGCTCGATTTGCAGCAGTCTATCACCCCAGTCGCCAAGAACCTCGCGCTGAGGTAGTGTAAACTAAGGGCTAGTGTCACTGGAATAGGACATACCCAAAGATAATTGAAGGCATATGGAGCTACTTGCCCGACAATTACCTTGTAGTTGAGGGTACATAGGAAGTCCAACCGATAACAACGCGCATGCGTGACCTTCCTGTCTTCCCAGCACAATGCCCACGCTACTGAGAGTCTTGATTTCTGGAATCAGTACAGGTTTTAACCCGTAGAAGGTATTTTAACCTAGTAAAGCGAAAATCGGGAGCAATAA')

      l = readlines('dataset_34405_4.txt')
      longestCommonSubseq(l[0], l[1]) >> testVs('CGTACATTTGTCCTTGGGGCTAGTCACAAATGCTGGGGTGGTTGTCGCGTGCAGTTTAGGTCGCCTACGCCGTATAAGGTAAATCGCATGGCCGGCAAGCTCCATTGCGCTGGGTGATTTGTATCGCCGGATTTTCAGTCTTCAGACCATGGGCAAGACCCGTGTGCTCATGTTGCAACCTTGCACTCCCTACTCCAGAGGGCGGCGTTTAAGAAGTAGGAGTAAAATCTAAAGAGTGGTGCGTTTAGTTCCGGTTCATAGCTTTGCGTGGTCCCTTCACGAACCTTATAATGCGGACGCCTCGGTTGAGTTGAGTCCGCTACAGAGAAAAAGGTGCCGGTTCCTGCGGGTTAGTGACTGTGGTAAGCGTAAGAACCTACCTTCAGTTTAGGTATCGGCAAGCTGTGCCAAAGTACCTAATTTCCCAGGGTTCGATGCCAAGACTAGAGTCTGTCCATCTCCATACCGAAGGTGGTAGATGGTACCATCCATTTGTATGAGCAATCCGGCTAGGCACCCGACCCTTTTCAGCCGAGTACACAGTTCAACCCAAACAGTGTGGCCAACGAACTTCAAAAGCTCG')

   # Longest Path in DAG
   def _5():
      wg = parseWtGraph('''0 1 7
0 2 4
2 3 2
1 4 1
3 4 3''')
      paths = topologicalSort(wg)
      print('paths',paths)
      d, p = longestWeightedPath(paths, '0')
      print(paths, 'd',d,'p', p, 'wg',wg)
      res = pathArrowStr(p, ' ')
      f"{d}\n{res}" >> testVs('9\n0 2 3 4')

#       wg = parseWGraph('''0->1:1
# 0->3:10
# 1->2:1
# 2->3:1''')
#       paths = topologicalSort(wg)
#       d, p = longestWeightedPath(paths, '0')
#       res = pathArrowStr(p, '->')
#       f"{d}\n{res}" >> testVs('10\n0->3')

   def _6():
      l = readlines('dataset_245_7.txt')
      wg = parseWtGraph('\n'.join(l[2:]))
      paths = topologicalSort(wg)
      d, p = longestWeightedPath(paths, l[0])
      res = pathArrowStr(p, '->')
      f"{d}\n{res}" >> testVs('172\n0->1->4->5->7->12->13->16->17->22->27->32->36->43->46->49')

if run(2):
   chap5._2()
if run(3):
   chap5._3()
if run(4):
   chap5._4()
if run(5):
   chap5._5()
if run(6):
   chap5._6()