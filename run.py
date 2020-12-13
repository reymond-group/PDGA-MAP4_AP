from pdga import PDGA
from pdga import get_map4_distancefn, get_map4_fingerprintfn, get_ap_distancefn, get_ap_fingerprintfn
import sys
import numpy as np


if __name__ == "__main__":

    if len(sys.argv) == 2:
        seed = int(sys.argv[1])
    elif  len(sys.argv) > 2:
        print("Too many arguments")
        sys.exit(-1)
    else:
        seed = np.random.randint(100000000)
        #print("Use seed:",seed)

    query_smiles = "C[C@H]1CCC[C@@H]2[C@H](CC(/C(C)=C/C3=CSC(C)=N3)OC(C[C@H](O)C(C)(C)C([C@H](C)[C@H]1O)=O)=O)O2"
    query_name = "EpothilonA"
    
    # To use MAP4 as fingerprint:
    # fingerprintfn = get_map4_fingerprintfn(return_strings=True, radius=1)
    # distancefn = get_map4_distancefn(use_string=True)

    # To use the RDKit Atom pair as fingerprint:
    fingerprintfn = get_ap_fingerprintfn()
    distancefn = get_ap_distancefn()
    
    # if you want to avoid the possibility of metylation of the amide bond set methyl to False
    ga = PDGA(pop_size=200, mut_rate=0.5, gen_gap=0.5, query=query_smiles, sim_treshold=0.6, \
              porpouse="cyclic", folder=f"{query_name}_{seed}", fingerprintfn=fingerprintfn, distancefn=distancefn, query_name=query_name,\
              is_peptide_sequence=False, methyl=True, verbose=False, seed=seed)
    
    
    # if you want to exclude specific building blocks (for a complete list refer to: https://github.com/reymond-group/PeptideDesignGA)
    # ga.sequence_rng.exclude_buildingblocks(['Hyp', 'Orn','bAla','Gaba','a5a','a6a','a7a','a8a','a9a'])
    
    
    
    # debug:
    ga.write_param()

    ga.set_time_limit('01:00:00')

    ga.run()
    



