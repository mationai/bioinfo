import timeit
import os
from sys import argv
from extensions import *

# def histogram[T](L:list[T]) -> dict[T, int]:
#    res = dict[T, int]()
#    for i in L:
#       res[i] = res[i]+1 if i in res else 1
#    return res

# def runningSums[T](L:list[T]) -> list[T]:
#    return [sum(L[:i]) for i in range(len(L)+1)]

# def mergeOrderedSeqs(seqs:list):
# # def mergeOrderedSeqs[T](seqs:list[T]) -> T:
#    """ Concat adjacent seqs. Eg. abc, bcd -> abcd
#    """
#    if len(seqs) == 1:
#       return seqs[0]
#    if len(seqs) < 1:
#       print('WARN: mergeOrderedSeqs -> "" due to empty input seqs')
#       return T('')
#    res = seqs[0]
#    # n = len(res)-1
#    for i, dna in enumerate(seqs[1:]): 
#       if res[i+1:] == dna[:-1]:
#          res += dna[-1]
#    return res

def printRuntime(t0:float, deci=1):
   print('Elapsed time:',
      round(timeit.default_timer() - t0, deci))

def sepWithSp(L:list) -> str:
   return ' '.join([str(x) for x in L])

def sepWithNL(L:list) -> str:
   return '\n'.join([str(x) for x in L])

def run(n:int) -> bool:
   x = argv[1] if len(argv) > 1 else ''
   if '..' in x:
      low, high = [int(s.strip()) for s in x.split('..')]
      if low <= n <= high:
         print('\n'+str(n)+':')
         return True
   else:
      if x=='' or int(x)==n:
         print('\n'+str(n)+':')
         return True
   return False

def testIn(B:list, a):
# def testIn[T](a:T, B:list[T]):
   if any(a==x for x in B):
      print("✓ result as expected")
   elif len(B)==1:
      print("  result:", a)
      print("expected:", B[0])
   else:
      print("a not in B")

# def testVerbose[T](a:T, b:T):
#    if a == b:
#       print("✓ result as expected")
#       return
#    iMismatch = Ints()
#    for i, c in enumerate(a):
#       if len(b) > i and c != b[i]: 
#          iMismatch.append(i)
#    if len(iMismatch) > 10:
#       print(f"{len(iMismatch)} or more mismatches, eg.:")
#    cnt = 0
#    for i in iMismatch:
#       if cnt < 10:
#          print(f"{i}: got {a[i]}, expected {b[i]}")
#          cnt += 1

def testVs(exp, res):
# def testVs[T](exp:T, res:T):
   if exp == res:
      print("✓ result as expected")
   else:
      print("  result:", res)
      print("expected:", exp)

def testColl(A:list, B:list):
# def testColl[T](A:List[T], B:List[T]):
   matchedX = set[int]()
   matchedY = set[int]()
   for x, a in enumerate(A): 
      for y, b in enumerate(B): 
         if a == b and x not in matchedX and y not in matchedY:
            matchedX.add(x)
            matchedY.add(y)
   if len(matchedX)==len(A) and len(matchedY)==len(B):
      print("✓ result as expected")
   else:
      print("  result:", B)
      print("expected:", A)

def parseInts(s:str, sep=' ') -> [int]:#Ints:
   return [int(x) for x in s.split(sep) if x.strip()]

# def parseIntsList(S:Strs, sep=' ') -> list[Ints]:
#    return [parseInts(s) for s in S]

def parseFloats(s:str, sep=' ') -> [float]:#Floats:
   return [float(x) for x in s.split(sep) if x.strip()]

# def parseToStrs(s:str, sep=' ') -> Strs:
#    return [x for x in s.split(sep) if x.strip()]

def parseSeqsStr(s:str, sep=' ') -> [str]:#Seqs:
   return [str(x) for x in s.split(sep) if x.strip()]

def parseSeqStrs(S:[str]) -> [str]:
# def parseSeqStrs(S:Strs) -> Seqs:
   return [str(x) for x in S if x.strip()]

def parseProfile(d) -> Profile:
   return {str(key): parseFloats(d[i]) for i, key in enumerate('ACGT')}

# def parseMat(st:str, lineSep='\n', sep=' ') -> IntMat:
#    return [[int(i) for i in line.split(sep)] for line in st.strip().split(lineSep)]

# def parseMats(st:str, sep='-\n') -> list[IntMat]:
#    return [parseMat(s) for s in st.split(sep)]

# def pairsStrToStrPairs(st:str, sep='\n', pairSep='|') -> StrPairs:
#    res = StrPairs()
#    for s in st.split(sep):
#       a, b = s.split(pairSep)
#       res.append((a, b))
#    return res

# def pairStrsToStrPairs(strs:Strs, sep='\n', pairSep='|') -> StrPairs:
#    res = StrPairs()
#    for s in strs:
#       a, b = s.split(pairSep)
#       res.append((a, b))
#    return res

# def listToStr[T](L:list[T], sep=' ') -> str:
#    return sep.join([str(x) for x in L])

# def toDashed[T](L:list[T]) -> str:
#    return '-'.join([str(i) for i in L])
   
# def spaceDashed[T](LL:list[list[T]], sort=True) -> str:
#    """ Returns eg. [[1, 2], [3, 4]] to '1-2 3-4'
#    """
#    res = Strs()
#    for L in LL:
#       res.append('-'.join([str(i) for i in L]))
#    return ' '.join(sorted(res) if sort else res)

# def toPairsStrs(pairs:Pairs, sep='|') -> Strs:
#    """ Returns list of "(a1|b1)" representation of input pairs
#    """
#    return [f"({a}{sep}{b})" for a, b in pairs]

# def toPairsStr(pairs:Pairs, sort=True) -> str:
#    """ Returns sorted(default) "(a1|b1) (a2|b2) ..." representation of input pairs
#    """
#    strs = pairs |> toPairsStrs
#    if sort:
#       strs.sort()
#    return strs |> sepWithSp

def adjGraphStr(sg:Graph) -> str:
   res = [f"{x}: {' '.join(sg[x])}" for x in sg]
   return '\n'.join(res)

# def pathArrowStr[T](S:list[T], arrow=' -> ') -> str:
#    return arrow.join([str(s) for s in S])

# def pathsArrowStr[T](SS:list[list[T]], arrow=' -> ') -> str:
#    return '\n'.join([pathArrowStr(S, arrow) for S in SS])

# def rmSpaces(s:str) -> str:
#    return s.replace(' ', '')

def dPath(p:str):
   dir = os.path.dirname(__file__)
   return os.path.join(dir, '../data', p)

# def file(p:str) -> File:
#    return open(dPath(p), 'r')

# def read(p:str, sz=100000) -> str:
#    return open(dPath(p), 'r').read(sz)

def readlines(p:str) -> list[str]:
   f = open(dPath(p), 'r')
   return f.read().splitlines()

def product(nums:list):
# def product[T](nums:list[T]) -> T:
   res = nums[0]
   for n in nums[1:]:
      if n == 0:
         return 0
      res *= n
   return res

# def maxIdx[T](M:list[T], minVal:T) -> int:
#    iLen = len(M)
#    score = minVal
#    iMax = 0
#    for i in range(iLen):
#       if M[i] > score:
#          score = M[i]
#          iMax = i
#    return iMax

# def maxXY[T](M:list[list[T]], minVal:T) -> tuple[int, int]:
#    jLen = len(M)
#    iLen = len(M[0])
#    score = minVal
#    iMax, jMax = 0, 0
#    for j in range(jLen):
#       for i in range(iLen):
#          if M[j][i] > score:
#             score = M[j][i]
#             iMax, jMax = i, j
#    return iMax, jMax

# def maxIn[T](M:list[list[T]], minVal:T) -> tuple[int, int]:
#    i, j = maxXY(M, minVal)
#    return M[j][i]

# def zeros(xLen:int, yLen:int) -> IntMat:
#    return [[0 for i in range(xLen)] for j in range(yLen)]

# def ones(xLen:int, yLen:int, v=1) -> IntMat:
#    return [[v for i in range(xLen)] for j in range(yLen)]

# ## ----- lambda key functions -----
# def getSI1[T](x:tuple[str, int]) -> int:
#    return x[1]
