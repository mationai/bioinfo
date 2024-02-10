from extensions import *
from seq2graph import *
from helpers import *


def parseGraph(st:str, sep='\n', arrow=' -> ') -> Graph:
   res = Graph()
   for s in [s.strip() for s in st.split(sep)]:
      lf, rt = s.split(arrow)
      res[lf] = rt.split(',') if ',' in rt else [rt]
   return res

def parseWGraph(st:str, sep='\n', arrow='->', wtSep=':') -> WGraph:
   """ parse weight graph string to WGraph, Eg.
   0->1:1  0->3:10  1->2:1  2->3:1 =>
   0: [(1,1), (3,10)],  1:[(2,1)],  2:[(3,1)]
   """
   res = WGraph()
   for s in [s.strip() for s in st.split(sep)]:
      src, rt = s.split(arrow)
      dst, wt = rt.split(wtSep)
      res[src] = res[src]+[(dst, int(wt))] if src in res else [(dst, int(wt))]
   return res

def readGraphStr(p:str) -> Graph:
   s = '\n'.join(open(dPath(p), 'r').readlines())
   return parseGraph(s)

def copyWG(wg:WGraph) -> WGraph:
   return {k: v[:] for k, v in wg.items()}


def eulerianPath(sg:Graph, kind:str, start:str) -> Strs:
   """ See http://www.graph-magics.com/articles/euler.php
   kind one of:
   'cycle'
   'path'
   'universal'
   """
   stack = Strs()
   res = Strs()
   keys = list(sg.keys())
   nVertices = len(keys)
   paths = {key: vals[:] for key, vals in sg.items()}
   inCnt = {key: 0 for key in keys}
   outCnt = {key: 0 for key in keys}
   cur = ''

   # Set in/out counts map
   for key, vals in sg.items():
      outCnt[key] += len(vals)
      for val in vals:
         inCnt[val] = 1 if val not in inCnt else inCnt[val] + 1

   eqIOCnts = sum(1 for k in keys if inCnt[k]==outCnt[k])
   notEqIOKeys = [k for k in keys if inCnt[k] != outCnt[k]]
   inGTOutKeys = [k for k in notEqIOKeys if inCnt[k] > outCnt[k]]
   outGTInKeys = [k for k in notEqIOKeys if outCnt[k] > inCnt[k]]
   maxResLen = sum(len(vals) for key, vals in sg.items()) + 1

   if len(keys) == 0:
      return res
   elif len(keys) == 1:
      cur = keys[0]
   elif start != '':
      cur = start

   # 1. Choose starting vertex
   # If all vertices have same fanout as fanin
   elif eqIOCnts == nVertices:
      # choose one with more than 1 fanout
      for i, key in enumerate(keys):
         if len(sg[key]) > 1:
            cur = key
            break
      # or one with more than 1 fanin if exists, else choose any, eg. 1st one
      if cur == '':
         cur = inGTOutKeys[0] if len(inGTOutKeys) > 0 else keys[0]

   # If all but 2 vertices have same out-degree as in-degree
   elif eqIOCnts == nVertices - 2:
      if len(outGTInKeys)==len(inGTOutKeys):
         outGTkey = outGTInKeys[0]
         inGTkey = inGTOutKeys[0]
         # If 1 of the 2 vertices has fanout 1 > fanin and the other is opposite then
         #  set cur to vertex w/ fanout 1 > fanin 
         if inCnt[inGTkey]+1 == outCnt[inGTkey] and outCnt[outGTkey]+1 == inCnt[outGTkey]:
            cur = outGTkey

   # Else no euler cycle exists if cur isn't set
   if kind=='cycle' and cur == '':
      return res

   # Path finding, set cur to key that has fanout > fanin
   if cur == '' and kind != 'cycle':
      cur = outGTInKeys[0] if len(outGTInKeys) > 0 else keys[0]

   # 2. Repeat until current vertex has no more out-going edges and stack is empty:
   while len(res) < maxResLen and cur != '':
      nexts = Strs() if kind != 'cycle' and cur not in paths else paths[cur]
      # nexts = Seqs() if kind != 'cycle' and cur not in paths else paths[cur]

      # If current vertex has no out-edges, add it to cycle and
      # remove the last vertex from the stack and set it as the current one
      if len(nexts) > 0:
         stack.append(cur)
         cur = nexts.pop(0)
      # Else, add the vertex to the stack, take any of its out-vertices and
      # remove that out-edge and set that out-vertex as the current vertex
      else:
         res.append(cur)
         cur = stack.pop() if len(stack) > 0 else ''

   if cur != '':
      res.append(cur)
   while len(stack) > 0:
      res.append(stack.pop())
   res.reverse()
   return res

def mergeOrderedSeqs[T](seqs:list[T]) -> T:
   """ Concat adjacent seqs. Eg. abc, bcd -> abcd
   """
   if len(seqs) == 1:
      return seqs[0]
   if len(seqs) < 1:
      print('WARN: mergeOrderedSeqs -> "" due to empty input seqs')
      return T('')
   res = seqs[0]
   # n = len(res)-1
   for i, dna in enumerate(seqs[1:]): 
      if res[i+1:] == dna[:-1]:
         res += dna[-1]
   return res

def genomeFromSeqs(seqs:Seqs, kind:str, start:str) -> str: 
   """ Reconstruct genome from seqs with method and kind (see eulerianPath())
   """
   sg = Seq2Graph(seqs, 'kmers', 0).sg
   path = eulerianPath(sg, kind, start)
   dna = path >> mergeOrderedSeqs
   if kind != 'universal':
      return dna 
   k = len(seqs[0]) - 1
   i = k // 2
   isOdd = k % 2 == 1 
   res = dna[i:-i]
   return res[1:] if isOdd else res

def pairedComposition(ptn:str, k:int, d:int) -> Pairs:
   """ Returns pairs of ptn formed by its (k, d)-mer
   """
   s = list(slide(ptn, k))
   l = len(ptn) -k -d -2
   return [(x, y) for x,y in zip(s[:l], s[k+d:])]

def unzip(strPairs:StrPairs) -> tuple[Strs, Strs]: 
   """ a, b = zip(*sps) (doesn't work in str)
   """
   A = [sp[0] for sp in strPairs]
   B = [sp[1] for sp in strPairs]
   return A, B

def unzipToSeqs(strPairs:StrPairs) -> tuple[Seqs, Seqs]: 
   """ a, b = zip(*sps) (doesn't work in str)
   """
   A = [str(sp[0]) for sp in strPairs]
   B = [str(sp[1]) for sp in strPairs]
   return A, B

def pathEdgesFromNodes[T](S:list[T]) -> list[T]:
   """ Eg. [a, b, c] -> [ab, bc]
   """
   if len(S) < 2:
      return S
   return [s+S[i+1] for i, s in enumerate(S[:-1])]

def isTunnelNode(src:str, sg:Graph) -> bool:
   """ Return if src in sg has only 1 other node pointing to it
   and if ag[src] has only 1 node
   Q: Still a tunnel of another node points to src and other nodes?
   """
   if src not in sg or len(sg[src]) != 1:
      return False
   cnt = 0
   for vals in sg.values():
      cnt += sum([1 for v in vals if v == src])
   return cnt == 1

def maxNonBranchingPaths(sg:Graph) -> list[Strs]:
   """ Returns non-branching paths in sg
   """
   res = list[Strs]()
   seen = StrsSet()
   for src, vals in sg.items():
      if not isTunnelNode(src, sg):
         for val in vals:
            path = [src, val]
            v = val
            isTunnel = isTunnelNode(v, sg)
            while isTunnel:
               v = sg[v][0] if len(sg[v]) > 0 else ''
               if v != '':
                  path.append(v)
                  seen.add(v)
               isTunnel = isTunnelNode(v, sg)
            res.append(path)
            seen.add(val)
         seen.add(src)
   for src, vals in sg.items():
      if not src in seen:
         path = [src]
         v = vals[0] if len(vals) > 0 else ''
         while v != '' and v != src:
            path.append(v)
            seen.add(v)
            v = sg[v][0] if len(sg[v]) > 0 else ''
         if v == src:
            path.append(src)
         res.append(path)
         seen.add(src)
   return res

def contigsFromSeq(ptn:str, k:int) -> Strs:
   """ Generate Contig - long, segments of non-branching genome from a str 
   Eg. (see L2-2 1.4.4), k=3 (kmer or edge len)
   ptn = TAATGCCATGGGATGTT
   ->
   [TAAT TGTT TGCCAT ATG ATG ATG TGG GGG GGAT]
   """
   kmers = list(slide(ptn, k))
   sg = Seq2Graph(kmers, 'kmers', 0).sg
   path = eulerianPath(sg, 'path', '')
   stack = Strs()
   res = Strs()

   for p in path:
      if p in sg and len(sg[p]) == 1:
         stack.append(p)
      else:
         res.append(mergeOrderedSeqs(stack + [p]))
         stack = [p]
   return res

def contigs(kmers:Seqs) -> Strs:
   """ Generate Contig - long, segments of non-branching genome from kmers
   """
   sg = Seq2Graph(kmers, 'kmers', 0).sg
   paths = maxNonBranchingPaths(sg)
   return [mergeOrderedSeqs(p) for p in paths]

def genomeFromOrderedPairs(strPairs:StrPairs, k:int, d:int) -> str:
   """ Construct genome from ordered pairs, aka "string spelled by gapped patterns"
   """
   A, B = unzip(strPairs) 
   ptns1 = strPairs[0][0] + ''.join([a[-1] for a in A[1:]])
   ptns2 = strPairs[0][1] + ''.join([b[-1] for b in B[1:]])
   l = len(ptns1)
   if ptns1[k+d:] != ptns2[:l-k-d]:
      print('WARN: No string spelled by the gapped patterns')
   return str(ptns1 + ptns2[-(k+d):])

def genomeFromPairs(strPairs:StrPairs, k:int, d:int) -> str:
   """ Construct genome from ordered pairs, aka "string spelled by gapped patterns"
   """
   pairsCnt = len(strPairs)
   A, B = unzipToSeqs(strPairs)
   a = genomeFromSeqs(A, 'path', '')
   target = a[k+d:]
   bStart = target[:k-1]
   b = genomeFromSeqs(B, 'path', bStart)
   tried = StrsSet()
   aStart = a[:k-1]
   resLen = pairsCnt + k + d + k-1
   iA = 0 

   while len(tried) < pairsCnt:
      for i in range(len(a)):
         b_ = b[i:]+b[:i]
         if target == b_[:-k-d]:
            res = a + b_[-k-d:]
            if len(res) == resLen:
               return a + b_[-k-d:]

      res = a + b[-k-d:]
      if len(res) == resLen:
         return a + b[-k-d:]
      tried.add(aStart)
      aStart = str(A[iA][:k-1])
      a = genomeFromSeqs(A, 'path', aStart)
      target = a[k+d:]
      bStart = target[:k-1]
      b = genomeFromSeqs(B, 'path', bStart)
      iA += 1 
   return ''
