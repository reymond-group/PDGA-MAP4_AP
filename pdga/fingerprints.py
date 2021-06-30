
from map4 import MAP4Calculator
import numpy as np
from rdkit.Chem import DataStructs
from rdkit.Chem.AtomPairs import Pairs

def tanimoto_map4_string(a, b):
    a = set(a)
    b = set(b)
    return 1.0 - (len(a.intersection(b))/len(a.union(b)))

def dice_map4_string(a, b):
    a = set(a)
    b = set(b)
    intersect = len(a.intersection(b))
    return 1 - ( 2 * intersect )/ (len(a) + len(b))

def minhash_map4_distance(a, b):
    return 1.0 - np.float(np.count_nonzero(a == b)) / np.float(len(a))


class get_map4_fingerprintfn:
    def __init__(self,dimensions=1024, return_strings=True, radius=1):
        self.MAP4 = MAP4Calculator(dimensions=dimensions, return_strings=return_strings, radius=radius)
        self.dimensions = dimensions
        self.return_strings = return_strings
        self.radius = radius

    def __call__(self, *args, **kwargs):
        return self.MAP4.calculate_many(*args, **kwargs)

    def __repr__(self) -> str:
        if self.return_strings:
            return f"MAP4 Unhashed: get_map4_fingerprintfn(radius={self.radius}, return_strings={self.return_strings})"
        else:
            return f"MAP4 MinHashed: get_map4_fingerprintfn(dimensions={self.dimensions}, radius={self.radius})"

class get_map4_distancefn:
    def __init__(self, use_string=True, tanimoto=True):
        self.use_string = use_string
        self.tanimoto = tanimoto
        if not use_string:
            self._fn = minhash_map4_distance
        else:
            if tanimoto:
                self._fn = tanimoto_map4_string
            else:
                self._fn = dice_map4_string

    def __call__(self, a, b):
        return self._fn(a,b)        

    def __repr__(self) -> str:
        if self.use_string:
            if self.tanimoto:
                return f"MAP4 Unhashed Tanimoto Distance: get_map4_distancefn(use_string={self.use_string}, tanimoto={self.tanimoto})"
            else:
                return f"MAP4 Unhashed Dice Distance: get_map4_distancefn(use_string={self.use_string}, tanimoto={self.tanimoto})"
        else:
            return f"MAP4 Minhashed Estimated Tanimoto Distance: get_map4_distancefn(use_string={self.use_string})"


class get_ap_fingerprintfn:
    def __call__(self, mols):
        return [Pairs.GetAtomPairFingerprint(x) for x in mols]

    def __repr__(self) -> str:
        return f"AtomPairFingerprint: get_ap_fingerprintfn()"


class get_ap_distancefn:
    def __call__(self, a, b):
        return 1 - DataStructs.DiceSimilarity(a,b)

    def __repr__(self) -> str:
        return f"AtomPairFingerprint Dice Distance: get_ap_distancefn()"
