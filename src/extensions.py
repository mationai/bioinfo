from typing import List, Set, Tuple, Dict

Ints = list#List[int]
Seqs = list#List[str]
Strs = list#List[str]
Floats = list#List[float]
SeqsSet = Set[str]
StrsSet = Set[str]
GraphDict = Dict[str, str]
Profile = Dict[str, List[float]]
SeqFloatDict = Dict[str, float]
SeqIntDict = Dict[str, int]
IntIntDict = Dict[int, int]
Graph = dict#Dict[str, Strs]
WtGraph = Dict[str, List[Tuple[str,int]]]
WPaths = List[Tuple[str, List[Tuple[str,int]]]]
Pairs = List[Tuple[str, str]]
SIPair = Tuple[str, int]
StrPair = Tuple[str, str]
StrPairs = List[StrPair]
IntMat = List[Ints]
StrMat = List[Strs]
StrMatMap = Dict[str, StrMat]
IntMatMap = Dict[str, IntMat]
ScoreMat = Dict[str, Dict[str,int]]

# extend str:
#    def __add__(s1:str, s2:str) -> str:
#       return str(str(s1) + str(s2))

#    def reverse(s1:str) -> str:
#       st = str(s1)
#       return str(''.join([s for s in reversed(st)]))

#    def replace(s1:str, frm:str, to:str) -> str:
#       st = str(s1)
#       return str(st.replace(frm, to))

#    def __sub__(s1:str, s2:str) -> int:
#       smLen = len(s1) if len(s1) < len(s2) else len(s2)
#       diff = len([1 for i in range(smLen) if s1[i] != s2[i]])
#       return diff + abs(len(s1) - len(s2)) 

# extend List[T]:
#    def sortedByList(X, Y):
#       return [x for _,x in sorted(zip(Y,X))]

# extend str:
#    # eslice additions: https://github.com/str-lang/str/commit/ae1384cdab890f888fe751d515bc50fd385726a8
#    def _make_from_range(self: str, start: int, stop: int, step: int, length: int):
#       p = ptr[byte](length)
#       j = 0
#       for i in range(start, stop, step):
#          p[j] = self.ptr[i]
#          j += 1
#       return str(p, length)

#    def __getitem__(self: str, s: esslice):
#       start, stop, step, length = slice.adjust_indices(len(self), step=s.step)
#       return self._make_from_range(start, stop, step, length)

#    def __getitem__(self: str, s: rsslice):
#       start, stop, step, length = slice.adjust_indices(len(self), start=s.start, step=s.step)
#       return self._make_from_range(start, stop, step, length)

#    def __getitem__(self: str, s: lsslice):
#       start, stop, step, length = slice.adjust_indices(len(self), end=s.end, step=s.step)
#       return self._make_from_range(start, stop, step, length)

#    def __getitem__(self: str, s: sslice):
#       start, stop, step, length = slice.adjust_indices(len(self), start=s.start, stop=s.end, step=s.step)
#       return self._make_from_range(start, stop, step, length)

   # end eslice additions

# extend dict:
#    def sortedKeysByVal(D):
#       T = sorted([v for v in D.values()])
#       return [x for _,x in sorted(zip(Y,X))]