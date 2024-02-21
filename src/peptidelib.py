from superpipe import pipes
from extensions import *
from helpers import *
from seqlib import *

CodonLU = {
   "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
   "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
   "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
   "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
   "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
   "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
   "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
   "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
   "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
   "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
   "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
   "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
   "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
   "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
   "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
   "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",
}
CodonRevLU = {
   # Ala
   'A': 'GCU GCC GCA GCG'.split(' '),
   # Cys
   'C': 'UGU UGC'.split(' '),
   # Asp
   'D': 'GAU GAC'.split(' '),
   # Glu
   'E': 'GAA GAG'.split(' '),
   # Phe
   'F': 'UUU UUC'.split(' '),
   # Gly
   'G': 'GGU GGC GGA GGG'.split(' '),
   # His
   'H': 'CAU CAC'.split(' '),
   # Ile
   'I': 'AUU AUC AUA'.split(' '),
   # Lys
   'K': 'AAA AAG'.split(' '),
   # Leu
   'L': 'CUU CUC CUA CUG UUA UUG'.split(' '),
   # Met
   'M': ['AUG'],
   # Asn
   'N': 'AAU AAC'.split(' '),
   # Pro
   'P': 'CCU CCC CCA CCG'.split(' '),
   # Gln
   'Q': 'CAA CAG'.split(' '),
   # Arg
   'R': 'CGU CGC CGA CGG AGA AGG'.split(' '),
   # Ser
   'S': 'UCU UCC UCA UCG AGU AGC'.split(' '),
   # Thr
   'T': 'ACU ACC ACA ACG'.split(' '),
   # Val
   'V': 'GUU GUC GUA GUG'.split(' '),
   # Trp
   'W': ['UGG'],
   # Tyr
   'Y': 'UAU UAC'.split(' '),
   '*': 'UAG UGA UAA'.split(' '),
}
AbrevLU = {
   'Tyr': 'Y',
   'Val': 'V',
   'Gly': 'G',
   'Ala': 'A',
   'Glu': 'E',
   'Asp': 'D',
   'Leu': 'L',
   'Arg': 'R',
   'Pro': 'P',
   'Gln': 'Q',
   'His': 'H',
   'Met': 'M',
   'Ile': 'I',
   'Arg': 'R',
   'Ser': 'S',
   'Thr': 'T',
   'Lys': 'K',
   'Asn': 'N',
   'Phe': 'F',
   'Trp': 'W',
   'Cys': 'C',
}
CodonKeyLen = 3


def abreviate(s:str, sep='-') -> str:
   return ''.join([AbrevLU[r] for r in s.split(sep)])

def strToRNA(dna:str) -> str:
   return dna.replace('T', 'U')
def seqToRNA(dna:str) -> str:
   return dna.replace('T', 'U')

def strToDNA(dna:str) -> str:
   return dna.replace('U', 'T')
def seqToDNA(dna:str) -> str:
   return dna.replace('U', 'T')

def translateProtein(ptn:str) -> str:
   """ Return mapping of ptn - str of CodonTable keys
   """
   if len(ptn) % 3 != 0:
      print("WARN: translateProtein() requires input of length % 3")
      return ''
   lu = CodonLU
   return ''.join([lu[str(c)] for c in slide(ptn, 3, 3) if lu[str(c)] != '*'])

def permCountsOf(ptn:str) -> int:
   """ Return count of possible keys mapping of ptn (str of CodonTable values)
   """
   return product([len(CodonRevLU[c]) for c in ptn])

def encodePeptide(dna:str, ptn:str) -> Strs:
   """ Return Peptide encoding for ptn in dna
   """
   res = Strs()
   l = CodonKeyLen
   ptnLen = len(ptn)*l
   for Ptn in slide(dna, ptnLen, 1):
      pVs = ''
      rVs = ''
      P = seqToRNA(Ptn)
      R = seqToRNA(reverseComp(Ptn))
      # print(Ptn)
      for p, r in [(P[i:i+l], R[i:i+l]) for i in range(0, ptnLen, l)]:
         if str(p) in CodonLU:
            pVs += CodonLU[str(p)]
         if str(r) in CodonLU:
            rVs += CodonLU[str(r)]
      # print(pVs, ptn, rVs)
      if pVs==ptn or rVs==ptn:
         res.append(str(P))
   return [strToDNA(s) for s in res]

def cycloSubpeptidesCount(l:int) -> int:
   """ Returns # of subpetitdes for length l, including cycles; excludes empty and itself
   """
   return l * (l - 1)

def linearSubpeptidesCount(l:int) -> int:
   """ Returns # of subpetitdes for length l, excluding cycles; includes empty and itself
   """
   return cycloSubpeptidesCount(l) +2 -sum(range(1, l-1)) # range() lists subpeptides counts that wraps

AminoMass = {
   'G': 57,
   'A': 71,
   'S': 87,
   'P': 97,
   'V': 99,
   'T': 101,
   'C': 103,
   'I': 113,
   'L': 113,
   'N': 114,
   'D': 115,
   'K': 128,
   'Q': 128,
   'E': 129,
   'M': 131,
   'H': 137,
   'F': 147,
   'R': 156,
   'Y': 163,
   'W': 186,
} # 20 amino acids
Aminos = list(AminoMass.keys())
AminoMasses = list(set(AminoMass.values()))
UniqueMassAminos = [a for a in Aminos if a not in ['L','Q']]

# "Extended" non-proteinogenic amino acids, expanding the number of possible
#  building blocks for antibiotic peptides from 20 to 144 (mass 57 - 200).
ExtAminoMass = {chr(i): i for i in range(57, 201)}
ExtAminoRevLU = {i: chr(i) for i in range(57, 201)}
ExtAminoMasses = list(ExtAminoMass.values())
ExtAminos = list(ExtAminoMass.keys())


def peptideMasses(peptide:str, massLU=AminoMass) -> Ints:
   """ Returns list of masses mapped from peptide strs
   """
   return [massLU[s] for s in peptide]

def peptideMass(peptide:str, massLU=AminoMass) -> int:
   """ Returns mass of sum of strings in peptide
   """
   return sum(peptideMasses(peptide, massLU))

@pipes
def linearSpectrum(peptide:str, massLU=AminoMass) -> Ints:
   """ Returns collection of all of the masses of its nonwrapping subpeptides, with
   masses ordered from smallest to largest. 
   """
   mass = peptideMasses(peptide, massLU) >> runningSums
   l = len(peptide)
   # self-referencing nested ("triangular") loop to iterate through what mass of j - i
   # , which should computes mass of linear subpeptides. 
   # cycloSpectrum() with this result + masses of wrapSubpeptides() returns few masses short
   #  however, and wrapSubpeptides() is verified correct, so perhaps this doesn't compute 
   #  all linear subpeptide masses.
   res = [0] + [mass[j] - mass[i] for i in range(l) for j in range(i+1, l+1)]
   return sorted(res)

def subpeptides(peptide:str) -> Strs:
   """ Returns all subpeptides
   """
   pLen = len(peptide)
   res = Strs()
   wrapped = peptide + peptide
   for start in range(0, pLen):
      for length in range(1, pLen):
         res.append((wrapped[start:start + length]))
   res.append(peptide)
   return res

def cycloSpectrum(peptide:str, massLU=AminoMass) -> Ints:
   """ Like linearSpectrum, but includes masses of wrapping subpeptides.
   """
   masses = [0] + [peptideMass(p, massLU) for p in subpeptides(peptide)]
   return sorted(masses)

def expand(peptides=[''], expansions=Aminos) -> Strs:
   """ Eg. peptides=[ab bc], expansions=[a d] -> [aba abd bca bcd]
   """
   return [p + c for p in peptides for c in expansions]

def isConsistent(peptide:str, spectrum:Ints, massLU=AminoMass) -> bool:
   """ A linear peptide is consistent with spectrum if all mass in its
    theoretical spectrum is contained in spectrum
   """
   S = spectrum[:]
   for i in linearSpectrum(peptide, massLU):
      if i not in S:
         return False
      S.remove(i)
   return True

def countPeptidesFrom(mass:int, masses=dict()) -> tuple: #[int, IntIntDict]:
   """ Returns count of all possible peptides whose mass is mass
   """
   if mass == 0:
      return 1, masses
   elif mass < 57:
      return 0, masses
   elif mass in masses:
      return masses[mass], masses

   res = 0
   for i in AminoMasses:
      n, masses = countPeptidesFrom(mass - i, masses)
      res += n
   masses[mass] = res
   return res, masses


def cyclopeptidesSequencing(spectrum:Ints) -> list: #[Ints]:
   """ Cyclopeptide Sequencing
   A brute force way is to check if spectrum == cycloSpectrum(peptide) for all peptide
    where mass of peptide == parentMass of spectrum. This will take a while however.

   Via Branch-and-Bound:
   If a mass appears more than once in the theoretical spectrum of the linear peptide,
    then it must appear at least that many times in Spectrum in order for the
    linear peptide to be consistent with Spectrum.

   The key to the algorithm is that every linear subpeptide of a cyclic peptide Peptide
    is consistent with Cyclospectrum(Peptide). Thus, to solve the Cyclopeptide Sequencing
    Problem for Spectrum, we can safely ban all peptides that are inconsistent with Spectrum
    from the growing set Peptides, which powers the bounding step that we described above.

   For example, the linear peptide VKF (with spectrum {0, 99, 128, 147, 227, 275, 374})
    will be banned because it is inconsistent with Tyrocidine B1’s spectrum.
    But the linear peptide VKY will not be banned because every mass in its theoretical
    spectrum ({0, 99, 128, 163, 227, 291, 390}) is present in Tyrocidine B1's spectrum.

   Note this method only works in ideal spectrums, ie. when the spectrum of a peptide
    coincides exactly with its theorectical spectrum.
   """
   expansions = [p for p in expand() if isConsistent(p, spectrum)]
   # print('expansions', expansions)
   parentMass = spectrum[-1]
   peptides = set([''])
   res = set()
   toDel = set()

   while len(peptides) > 0:
      peptides = set(expand(list(
         peptides.difference(toDel)), expansions))
      for peptide in peptides:
         if peptideMass(peptide) == parentMass:
            if cycloSpectrum(peptide) == spectrum:
               res.add('-'.join(
                  [str(p) for p in peptideMasses(peptide)]))
            toDel.add(peptide)
         elif not isConsistent(peptide, spectrum):
            toDel.add(peptide)
   return list(res)

def peptideScore(peptide:str, spectrum:Ints, kind='linear', massLU=AminoMass) -> int:
   """ Count (with min()) masses common in both spectrum and in cycloSpectrum of peptide.
   'linearScore' if kind is 'linear', cycloScore otherwise
   """
   linSpectrum = linearSpectrum(peptide, massLU) if kind=='linear' else cycloSpectrum(peptide, massLU)
   massCnt = histogram(linSpectrum)
   expMassCnt = histogram(spectrum)
   return sum([min(expMassCnt[m], massCnt[m]) for m in massCnt if m in expMassCnt])

def trimPeptides(peptides:Strs, spectrum:Ints, topN:int, massLU=AminoMass, minScore=0) -> Strs:
   """ Sorts all leaderboard peptides according to their scores.
   Return all peptides with score better or equal to score of topN
   If minScore is passed, it will be used instead of topN to obtain the min score.
   """
   if len(peptides) <= topN and minScore==0:
      return peptides
   scores = [peptideScore(p, spectrum, 'linear', massLU) for p in peptides]
   scoredRes = [(p, peptideScore(p, spectrum, 'linear', massLU)) for p in peptides]
   sortedRes = sorted(scoredRes, key=getSI1)[::-1]
   minScore = sortedRes[topN-1][1] if minScore==0 else minScore
   return [x[0] for x in sortedRes if x[1] >= minScore]

def topCyclopeptideSequencing(spectrum:Ints, topN:int) -> Ints:
   """ aka Leaderboard Cyclopeptide Sequencing
   Returns top peptide masses
   """
   parentMass = spectrum[-1]
   topPeptides = set([''])
   topPeptide = ''
   toDel = set()

   while len(topPeptides) > 0:
      leaderboard = set(expand(
         topPeptides.difference(toDel)))
      for p in leaderboard:
         mass = peptideMass(p)
         if mass == parentMass:
            if peptideScore(p, spectrum) > peptideScore(topPeptide, spectrum):
               topPeptide = p
         elif mass > parentMass:
            # leaderboard.remove(p)
            toDel.add(p)
      topPeptides = set(trimPeptides(list(leaderboard), spectrum, topN))
   return peptideMasses(topPeptide)[::-1]

def topPeptidesSequencing(spectrum:Ints, aminoKind='', topN=1000, altAminos=Strs()) -> list: #[Ints]:
   """ Like Leaderboard Cyclopeptide Sequencing, but
   Returns all top peptides masses
   """
   parentMass = spectrum[-1]
   topPeptides = ['']
   leaderboard = set(topPeptides)
   aminos = ExtAminos if aminoKind=='extended' else UniqueMassAminos
   massLU = ExtAminoMass if aminoKind=='extended' else AminoMass
   alphabet = altAminos if len(altAminos) > 0 else aminos

   while len(leaderboard) > 0:
      leaderboard = set(expand(leaderboard, alphabet))
      for p in leaderboard:
         mass = peptideMass(p, massLU)
         if mass == parentMass:
            score = peptideScore(p, spectrum, 'cyc', massLU)
            topScore = peptideScore(topPeptides[0], spectrum, 'cyc', massLU)
            if score >= topScore:
               topPeptides = [p] if score > topScore else topPeptides+[p]
               # print 'topScore', topScore, 'topPeps', topPeptides
         elif mass > parentMass:
            leaderboard.remove(p)
      leaderboard = set(trimPeptides(list(leaderboard), spectrum, topN, massLU))
   return [peptideMasses(p, massLU) for p in topPeptides]

def spectralConvolution(spectrum:Ints, gte=57, lte=200) -> Ints:
   """ 
   Returns list of elements in the convolution of Spectrum.
    If an element has multiplicity k, it should appear exactly k times;
   """
   return [i-j for i in spectrum for j in spectrum if gte <= i-j <= lte]

def getTopMasses(L:list, topN:int) -> Ints: # L:[SIPair], topN:int) -> Ints:
   if len(L) <= topN:
      return [int(x[0]) for x in L]
   minCnt = L[topN-1][1]
   return [int(x[0]) for x in L if x[1] >= minCnt]

@pipes
def convolutionSequencing(spectrum:Ints, topM=20, topN=60) -> list: #[Ints]:
   counts = spectralConvolution(spectrum) >> histogram
   massCnts = [(str(k), v) for k, v in counts.items()]
   sortedMasses = sorted(massCnts, key=getSI1)[::-1]
   topMasses = getTopMasses(sortedMasses, topM)
   alphabet = [ExtAminoRevLU[i] for i in topMasses]
   return topPeptidesSequencing(spectrum, 'extended', topN, alphabet)

def floatsSpectrumToInts(spectrum:Floats) -> Ints:
   """ spectrum are ions assumed to have a +1 charge and with a mass error of 0.3 Da.
    Outputs an integer mass spectrum.
   """
   masses = [x-1.0 for x in spectrum] # remove +1 charge
   rounded = Ints()
   for mass in masses:
      if mass - int(mass) == 0.5: # rounded up & down if x.5
         rounded.append(int(mass)+1)
         rounded.append(int(mass))
      else:
         # take 2 common ints from 3 of: -.3, +0, +.3 (all +.5 to floor)
         m3 = [int(mass+0.2), int(mass+0.5), int(mass+0.8)]
         mid = sorted(m3)[1]
         rounded.append(mid)
   return sorted(rounded)
    
@pipes
def protonMassSpecCyclopeptideSequencing(fSpectrum:Floats, topM:int, topN:int) -> Ints:
   """ Didn't produce correct answer
   """
   spectrum = floatsSpectrumToInts(fSpectrum)
   counts = spectralConvolution(spectrum) >> histogram
   massCnts = [(str(k), v) for k, v in counts.items()]
   sortedMasses = sorted(massCnts, key=getSI1)[::-1]
   topMasses = getTopMasses(sortedMasses, topM)
   alphabet = [ExtAminoRevLU[i] for i in topMasses]
   tops = topPeptidesSequencing(spectrum, 'extended', topN, alphabet)
   return tops[0]


def wrapSubpeptides(peptide:str) -> Strs:
   """ Returns only wrapping subpeptides
   Works, but not used
   """
   pLen = len(peptide)
   if pLen < 2:
      return Strs()
   res = set() #StrsSet()

   for i in range(1, len(peptide)-2):
      res.add(peptide[-1]+peptide[:i+1])
      res.add(peptide[-i-1:]+peptide[0])
   res.add(peptide[-1]+peptide[0])
   return list(res)
