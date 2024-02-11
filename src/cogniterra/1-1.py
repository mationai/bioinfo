from superpipe import pipes
import sys
import timeit
sys.path.append('src')
from helpers import *
from seqlib import *

@pipes
class chap1:
   # Get all k-mers in seq (including repeats)
   def _1():
      slide('CAATCCAAC', 5) >> sepWithSp\
         >> testVs('CAATC AATCC ATCCA TCCAA CCAAC')

   # Find occurences count of ptn in dna
   def _2():
      indexesIn('GCGCG', 'GCG') >> len >> testVs(2)

      [dna, ptn] = readlines('dataset_33716_2a.txt')
      indexesIn(dna, ptn) >> len >> testVs(26)

      [dna, ptn] = readlines('dataset_33716_2.txt')
      indexesIn(dna, ptn) >> len >> testVs(25)

   # Find most freq k-mers in dna
   def _3():
      [kmrs, cnt] = mostFreqKmers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4)
      print('cnt:', cnt)
      kmrs >> sepWithSp >> testVs('CATG GCAT')

      [txt, cnt] = readlines('dataset_33716_3a.txt')
      [kmrs, cnt] = mostFreqKmers(txt, int(cnt))
      print('cnt:', cnt)
      kmrs >> sepWithSp >> testVs('ATATTCCACCCTTG ATTCCACCCTTGAT TATTCCACCCTTGA TTCCACCCTTGATT')

      [txt, cnt] = readlines('dataset_33716_3.txt')
      [kmrs, cnt] = mostFreqKmers(txt, int(cnt))
      print('cnt:', cnt)
      kmrs >> sepWithSp >> testVs('GCCTATCCGATAAG')

   # Find reverse complement of dna
   def _4():
      rc = reverseComp('AAAACCCGGT')
      rc >> testVs('ACCGGGTTTT')

      dna = readlines('dataset_33716_4a.txt')[0]
      rc = reverseComp(dna)
      [rc[:6], rc[-6:]] >> testVs(['GGGTTC', 'CGATGT'])

      dna = readlines('dataset_33716_4.txt')[0]
      rc = reverseComp(dna)
      [rc[:6], rc[-6:]] >> testVs(['GTGGAC', 'GGATGG'])

   # Find Indices of ptn in dna
   def _5():
      indexesIn('GATATATGCATATACTT', 'ATAT') >> sepWithSp\
         >> testVs('1 3 9')

      [ptn, dna] = readlines('dataset_33716_5a.txt')
      indexesIn(dna, ptn) >> sepWithSp\
         >> testVs('22 41 81 93 124 186 223 280 305 388 395 404 411 418 443 498 507 536 600 659 715 778 785 839 846 853 935 978 1016 1048 1063 1082 1156 1171 1203 1210 1319 1326 1377 1384 1401 1435 1442 1449 1492 1564 1616 1644 1651 1658 1665 1698 1714 1738 1845 1925 1936 1998 2134 2155 2246 2331 2346 2445 2488 2524 2599 2638 2719 2726 2752 2780 2858 2874 2899 2937 2983 2990 3001 3024 3048 3092 3175 3203 3215 3282 3306 3332 3411 3418 3442 3449 3464 3528 3535 3552 3577 3606 3676 3704 3728 3786 3813 3820 3835 3842 3870 3916 3961 4003 4036 4069 4076 4101 4129 4136 4196 4203 4243 4315 4381 4388 4472 4498 4505 4543 4569 4626 4705 4751 4758 4789 4844 4851 4908 4923 4930 4937 5042 5049 5082 5149 5172 5219 5236 5252 5323 5330 5360 5442 5457 5464 5493 5500 5539 5562 5648 5674 5689 5822 5845 5852 5868 5928 5949 6008 6063 6070 6101 6108 6128 6179 6213 6228 6311 6318 6346 6441 6452 6511 6648 6759 6814 6829 6836 6843 6882 6889 6908 7003 7049 7056 7086 7122 7154 7161 7168 7187 7283 7301 7308 7331 7338 7364 7383 7425 7450 7469 7476 7522 7583 7615 7669 7676 7744 7812 7819 7826 7881 7888 7895 7926 8057 8107 8114 8129 8152 8167 8217 8233 8288 8323 8353 8419 8507 8540 8633 8652 8694 8739 8746 8756 8798 8805 8822 8839 8846 8853 8940 8961 9008 9015 9086 9111 9167 9176 9183 9199 9224 9253 9280 9293 9332 9348 9378 9438 9465 9472 9479 9530 9565 9594')

      [ptn, dna] = readlines('dataset_33716_5.txt')
      indexesIn(dna, ptn) >> sepWithSp\
         >> testVs('15 49 95 110 121 128 163 170 222 260 299 352 359 442 470 519 547 682 699 827 852 919 970 977 1000 1007 1107 1176 1209 1249 1261 1294 1313 1320 1400 1419 1540 1569 1585 1618 1649 1675 1682 1699 1719 1734 1810 1867 1874 1881 1933 1960 1995 2002 2137 2144 2202 2209 2244 2251 2292 2318 2347 2378 2473 2632 2694 2718 2725 2732 2787 2824 2931 2979 3013 3041 3083 3093 3174 3259 3279 3294 3311 3349 3356 3363 3494 3501 3512 3540 3669 3684 3691 3741 3776 3783 3790 3797 3835 3860 3867 3874 3903 3921 3930 3937 3975 4005 4012 4035 4138 4202 4302 4352 4440 4521 4592 4621 4652 4702 4709 4723 4730 4784 4791 4836 4865 4918 4948 4985 4992 5031 5052 5141 5148 5228 5268 5379 5397 5414 5462 5493 5583 5617 5678 5685 5695 5720 5727 5858 5885 5892 5953 6002 6025 6042 6098 6114 6140 6220 6238 6276 6283 6396 6412 6436 6443 6450 6471 6531 6538 6672 6716 6736 6772 6783 6897 6929 6948 6963 7025 7032 7087 7094 7159 7224 7231 7246 7264 7291 7298 7351 7358 7382 7407 7414 7457 7513 7520 7553 7596 7614 7621 7628 7694 7701 7708 7717 7749 7766 7821 7846 7862 7935 7950 8010 8032 8076 8110 8129 8136 8143 8265')
   
   # Find Indices of ptn in Vibrio_cholerae.txt
   def _5b():
      dna = readlines('Vibrio_cholerae.txt')[0]
      indexesIn(dna, 'CTTGATCAT') >> sepWithSp\
         >> testVs('60039 98409 129189 152283 152354 152411 163207 197028 200160 357976 376771 392723 532935 600085 622755 1065555')

   # Find Clumps - k-mers forming (L, t)-clumps in dna
   def _6():
      dna = 'CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA'
      getClumps(dna, 5, 50, 4) >> sepWithSp >> testVs('CGACA GAAGA')

      [dna, nums] = readlines('dataset_33716_6a.txt')
      [k, L, n] = parseInts(nums)
      getClumps(dna, k, L, n) >> sepWithSp\
         >> testVs('AAAGCGCAA AAGCGCAAA CAAAGCGCA')

      [dna, nums] = readlines('dataset_33716_6.txt')
      [k, L, n] = parseInts(nums)
      getClumps(dna, k, L, n) >> sepWithSp\
         >> testVs('AAATAGGGGG AATAGGGGGA ACTACCGAAT ACTCGCATGC CACACTAACG CCGGGTTTGT CCTGAGACGA GAAACTACGA GCCTCGCCTG GGAGATTACA TACCAAACCC TACTGATGGT')

   # Find length of clumps in E-coli
   def _6b():
      dna = readlines('E_coli.txt')[0]
      getClumps(dna, 9, 500, 3) >> len >> testVs(1904)

if run(1):
   chap1._1()
if run(2):
   chap1._2()
if run(3):
   chap1._3()
if run(4):
   chap1._4()
if run(5):
   chap1._5()
   t0 = timeit.default_timer()
   chap1._5b()
   printRuntime(t0)
   # mbp: ~5secs
if run(6):
   chap1._6()
   t0 = timeit.default_timer()
   chap1._6b()
   printRuntime(t0)
   # mbp: ~26secs