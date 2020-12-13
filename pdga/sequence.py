from numpy.core.fromnumeric import squeeze
from rdkit.Chem import rdmolfiles
from . import utils
import random 


B_SMILES = {'1': '[N:2][C@@H](C[N:2])[C:1](O)=O', '2': '[N:2][C@@H](CC[N:2])[C:1](O)=O',
            '3': '[N:2][C@@H](CCC[N:2])[C:1](O)=O', '4': '[N:2][C@@H](CCCC[N:2])[C:1](O)=O',
            '5': '[N:2][C@H](C[N:2])[C:1](O)=O', '6': '[N:2][C@H](CC[N:2])[C:1](O)=O',
            '7': '[N:2][C@H](CCC[N:2])[C:1](O)=O', '8': '[N:2][C@H](CCCC[N:2])[C:1](O)=O'}

AA_SMILES = {'A': '[N:2][C@@H](C)[C:1](O)=O', 'R': '[N:2][C@@H](CCCNC(N)=N)[C:1](O)=O',
                'N': '[N:2][C@@H](CC(N)=O)[C:1](O)=O', 'D': '[N:2][C@@H](CC(O)=O)[C:1](O)=O',
                'C': '[N:2][C@@H](CS)[C:1](O)=O', 'Q': '[N:2][C@@H](CCC(N)=O)[C:1](O)=O',
                'E': '[N:2][C@@H](CCC(O)=O)[C:1](O)=O', 'G': '[N:2]C[C:1](O)=O',
                'H': '[N:2][C@@H](CC1=CNC=N1)[C:1](O)=O', 'I': '[N:2][C@@H]([C@@H](C)CC)[C:1](O)=O',
                'K': '[N:2][C@@H](CCCCN)[C:1](O)=O', 'L': '[N:2][C@@H](CC(C)C)[C:1](O)=O',
                'M': '[N:2][C@@H](CCSC)[C:1](O)=O', 'F': '[N:2][C@@H](CC1=CC=CC=C1)[C:1](O)=O',
                'P': 'C1CC[N:2][C@@H]1[C:1](O)=O', 'S': '[N:2][C@@H](CO)[C:1](O)=O',
                'T': '[N:2][C@@H]([C@H](O)C)[C:1](O)=O', 'W': '[N:2][C@@H](CC1=CNC2=CC=CC=C12)[C:1](O)=O',
                'Y': '[N:2][C@@H](CC1=CC=C(C=C1)O)[C:1](O)=O', 'V': '[N:2][C@@H](C(C)C)[C:1](O)=O',
                'Ä': '[N:2][C@@H](C[S:1])[C:1](O)=O', 'Ö': '[N:2][C@@H](C[S:2])[C:1](O)=O',
                'Ü': '[N:2][C@@H](C[S:3])[C:1](O)=O',
                'Z': 'C1C(O)C[N:2][C@@H]1[C:1](O)=O',
                'O': '[N:2][C@@H](CCCN)[C:1](O)=O',
                'a': '[N:2][C@H](C)[C:1](O)=O', 'r': '[N:2][C@H](CCCNC(N)=N)[C:1](O)=O',
                'n': '[N:2][C@H](CC(N)=O)[C:1](O)=O', 'd': '[N:2][C@H](CC(O)=O)[C:1](O)=O',
                'c': '[N:2][C@H](CS)[C:1](O)=O', 'q': '[N:2][C@H](CCC(N)=O)[C:1](O)=O',
                'e': '[N:2][C@H](CCC(O)=O)[C:1](O)=O', 'g': '[N:2]C[C:1](O)=O',
                'h': '[N:2][C@H](CC1=CNC=N1)[C:1](O)=O', 'i': '[N:2][C@H]([C@@H](C)CC)[C:1](O)=O',
                'k': '[N:2][C@H](CCCCN)[C:1](O)=O', 'l': '[N:2][C@H](CC(C)C)[C:1](O)=O',
                'm': '[N:2][C@H](CCSC)[C:1](O)=O', 'f': '[N:2][C@H](CC1=CC=CC=C1)[C:1](O)=O',
                'p': 'C1CC[N:2][C@H]1[C:1](O)=O', 's': '[N:2][C@H](CO)[C:1](O)=O',
                't': '[N:2][C@H]([C@H](O)C)[C:1](O)=O', 'w': '[N:2][C@H](CC1=CNC2=CC=CC=C12)[C:1](O)=O',
                'y': '[N:2][C@H](CC1=CC=C(C=C1)O)[C:1](O)=O', 'v': '[N:2][C@H](C(C)C)[C:1](O)=O',
                'ä': '[N:2][C@H](C[S:1])[C:1](O)=O', 'ö': '[N:2][C@H](C[S:2])[C:1](O)=O',
                'ü': '[N:2][C@H](C[S:3])[C:1](O)=O',
                '!': '[N:2]CC[C:1](O)=O', '?': '[N:2]CCC[C:1](O)=O',
                '=': '[N:2]CCCC[C:1](O)=O', '%': '[N:2]CCCCC[C:1](O)=O',
                '$': '[N:2]CCCCCC[C:1](O)=O', '@': '[N:2]CCCCCCC[C:1](O)=O',
                '#': '[N:2]CC[C:1](O)=O'}


# list of possible aminoacids
AA = ['R', 'H', 'K', 'E', 'S', 'T', 'N', 'Q', 'G', 'P', 'A', 'V', 'I', 'L', 'F', 'Y', 'W', 'C', 'D', 'M', 'Z', 'O',
        '!', '?', '=', '%', '$', '@', '#']
# list of possible branching units (1=Dap, 2=Dab, 3=Orn, 4=Lys)
B = ['1', '2', '3', '4']
# list of possible C-terminals
CT = ['+']
# list of possible N-capping
NT = ['&']

interprete_dict = {'Arg': 'R', 'His': 'H', 'Lys': 'K', 'Asp': 'D', 'Glu': 'E', 'Ser': 'S', 'Thr': 'T', 'Asn': 'N',
                    'Gln': 'Q', 'Cys': 'C', 'Sec': 'U', 'Gly': 'G', 'Pro': 'P', 'Ala': 'A', 'Ile': 'I', 'Leu': 'L',
                    'Met': 'M', 'Phe': 'F', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Dap': '1', 'Dab': '2',
                    'BOrn': '3', 'BLys': '4', 'Hyp': 'Z', 'Orn': 'O', 'bAla': '!', 'Gaba': '?', 'dDap': '5',
                    'dDab': '6',
                    'dBOrn': '7', 'dBLys': '8', 'dArg': 'r', 'dHis': 'h', 'dLys': 'k', 'dAsp': 'd', 'dGlu': 'e',
                    'dSer': 's',
                    'dThr': 't', 'dAsn': 'n', 'dGln': 'q', 'dCys': 'c', 'dSec': 'u', 'dGly': 'g', 'dPro': 'p',
                    'dAla': 'a',
                    'dIle': 'i', 'dLeu': 'l', 'dMet': 'm', 'dPhe': 'f', 'dTrp': 'w', 'dTyr': 'y', 'dVal': 'v',
                    'dHyp': 'z', 'dOrn': 'o', 'a5a': '=', 'a6a': '%', 'a7a': '$', 'a8a': '@', 'a9a': '#',
                    'Cys1': 'Ä', 'Cys2': 'Ö', 'Cys3': 'Ü', 'dCys1': 'ä', 'dCys2': 'ö', 'dCys3': 'ü',
                    'Ac': '&', 'NH2': '+', 'met': '-', 'cy': 'X'}

interprete_rev_dict = {v: k for k, v in interprete_dict.items()}

T_SMILES = {'+': '[N:2]'}
C_SMILES = {'&': 'C[C:1](=O)'}

def is_cyclic(seq):
    if "X" in seq:
        return True
    else:
        return False

def add_cycle(seq):
    if not is_empty(seq) and not is_cyclic(seq):
        seq ='X' + seq
        for i in NT:
            seq = seq.replace(i, '')
        for i in CT:
            seq = seq.replace(i, '')    
    return seq


def remove_cycles(seq):
    seq = seq.replace("X", "")
    for i in NT:
        seq = seq.replace(i, '')
    for i in CT:
        seq = seq.replace(i, '')
    return seq


def ensure_maximum_one_cycle(seq):
    if is_cyclic(seq):
        seq = remove_cycles(seq)
        seq = add_cycle(seq)
    return seq


def linearize(seq):
    if is_cyclic(seq):
        return seq[1:]
    else:
        return seq

def is_empty(seq):
    if len(seq) == 0:
        return True
    if set(seq) <= set("X-"):
        return True
    else:
        False

def remove_constrained_activated_cystein_1(seq):
    seq = seq.replace("ÄÄ", '')
    return seq

def remove_constrained_activated_cystein_2(seq):
    seq = seq.replace("ÖÖ", '')
    return seq
    
def remove_constrained_activated_cystein_3(seq):
    seq = seq.replace("ÜÜ", '')
    return seq
    
def remove_constrained_activated_cystein(seq):
    seq = remove_constrained_activated_cystein_1(seq)
    seq = remove_constrained_activated_cystein_2(seq)
    seq = remove_constrained_activated_cystein_3(seq)
    return seq

def ensure_one_NT_at_beginning(seq):
    for i in NT:
        if not is_empty(seq):
            first_char = seq[0]
            if first_char == i:
                seq =  first_char + seq.replace(i,'')
            else:
                seq = seq.replace(i,'')
    return seq 

def ensure_one_CT_at_end(seq):
    for i in CT:
        if not is_empty(seq):
            last_char = seq[-1]
            if last_char == i:
                seq = seq.replace(i,'') + last_char
            else:
                seq = seq.replace(i,'')
    return seq 


def sanitize_C_and_N_Term(seq):
    if is_cyclic(seq):
        for i in NT:
            seq = seq.replace(i, '')
        for i in CT:
            seq = seq.replace(i, '')
    seq = ensure_one_CT_at_end(seq)
    seq = ensure_one_NT_at_beginning(seq)
    return seq


def sanitize_methylation(seq):
    double_methyl_present = "--" in seq
    while double_methyl_present:
        seq = seq.replace('--', '-')
        double_methyl_present = "--" in seq
    seq = remove_methylation_before_proline(seq)
    seq = remove_methylation_at_end_and_begining(seq)
    return seq



def sanitize_sequence(seq):
    seq = remove_constrained_activated_cystein(seq)
    seq = sanitize_C_and_N_Term(seq)
    seq = sanitize_methylation(seq)
    seq = ensure_maximum_one_cycle(seq)
    if is_empty(seq):
        return ""
    return seq 


def interprete(seq):
    """translates from 3letters code to one symbol
    
    Arguments:
        seq {string} -- 3 letters code seq (e.g. Ala-Gly-Leu)
    
    Returns:
        string -- one letter symbol seq (e.g. AGL)
    """

    new_seq = ''
    seq = seq.split('-')
    for bb in seq:
        new_seq += interprete_dict[bb]
    seq = new_seq
    return seq

def reinterprete(seq):
    """translates one symbol to three letters code
    
    Arguments:
        seq {string} -- one letter symbol seq (e.g. AGL)
    
    Returns:
        string -- 3 letters code seq (e.g. Ala-Gly-Leu)
    """

    new_seq = []
    for bb in seq:
        new_seq.append(interprete_rev_dict[bb])
    seq = '-'.join(new_seq)

    return seq



def split_seq_components(seq):
    """split seq in generations and branching units

    Arguments:
        seq {string} -- dendrimer sequence

    Returns:
        lists -- generations(gs, from 0 to..), branching units, terminal and capping
    """

    g = []
    gs = []
    bs = []
    t = []
    c = []

    for ix, i in enumerate(seq):
        if i not in ['1', '2', '3', '4', '5', '6', '7', '8']:
            if i in CT:
                t.append(i)
            elif i in NT:
                c.append(i)
            elif i == 'X':
                continue
            elif i == '-':
                if seq[ix - 1] in ['1', '2', '3', '4', '5', '6', '7', '8']:
                    bs.append(i)
                else:
                    g.append(i)
            else:
                g.append(i)
        else:
            gs.append(g[::-1])
            bs.append(i)
            g = []
    gs.append(g[::-1])
    gs = gs[::-1]
    bs = bs[::-1]

    return gs, bs, t, c




def find_aa_b_pos(seq):
    """finds aminoacids and branching unit positions in a given sequence
    
    Arguments:
        seq {string} -- peptide dendrimer sequence
    
    Returns:
        lists -- aminoacids and branching units positions, all position, terminal pos, capping 
    """

    aa = []
    b = []
    all_pos = []
    met = []

    for i, symbol in enumerate(seq):
        if symbol in ['X', 'Ö', 'Ü', 'Ä', 'ö', 'ü', 'ä'] + NT + CT:
            continue
        if symbol == '-':
            met.append(i)
            continue
        if symbol in B_SMILES.keys():
            b.append(i)
        elif symbol in AA_SMILES.keys():
            aa.append(i)
        all_pos.append(i)

    return aa, b, met, all_pos


def pick_aa_b_pos(seq, type_pos, verbose=False):
        """If type is aa, it returns an aminoacid position in the given sequence.
        if type is b, it returns a branching unit position in the given sequence.
        if type is all it returns a random position in the given sequence.
        
        Arguments:
            seq {string} -- peptide dendirmer sequence
            type_pos {string} -- aa, b or None
        
        Returns:
            int -- position
        """

        aa_pos, b_pos, met_pos, all_pos = find_aa_b_pos(seq)
        try:
            if type_pos == 'aa':
                return random.choice(aa_pos)
            elif type_pos == 'b':
                return random.choice(b_pos)
            elif type_pos == 'met':
                return random.choice(met_pos)
            elif type_pos == 'all':
                return random.choice(all_pos)
            else:
                if verbose:
                    print('not valid type, type has to be "aa", "b" or "all"')
        except:
            if verbose:
                print("Can't pick aa", seq)


def find_aminoacid_in_seq(seq):
    return pick_aa_b_pos(seq, 'aa')

def find_branching_in_seq(seq):
    return pick_aa_b_pos(seq, 'b')

def find_anything_in_seq(seq):
    return pick_aa_b_pos(seq, 'all')

def find_methylation_in_seq(seq):
    return pick_aa_b_pos(seq, 'met')



def form_SS(seq):
    """insertion of two ativated cys
    
    Arguments:
        seq {string} -- peptide seq
    
    Returns:
        string -- S-S cyclized peptide seq
    """

    act_cys = 'Ä'
    if 'Ä' in seq:
        act_cys = 'Ö'
        if 'Ö' in seq:
            act_cys = 'Ü'
            if 'Ü' in seq:
                return seq

    if len(seq.replace('X', '').replace('-', '')) <= 2:
        return seq

    # first active cys
    pos = find_aminoacid_in_seq(seq)
    seq_tmp = seq[:pos] + act_cys + seq[pos:]

    # second active cys
    pos = find_aminoacid_in_seq(seq)
    new_seq = seq_tmp[:pos] + act_cys + seq_tmp[pos:]

    # prevents to activated cys next to each other
    if act_cys + act_cys not in new_seq:
        seq = new_seq

    return seq


def methylate(seq):
    pos = find_anything_in_seq(seq)
    if not pos:
        return seq
    if pos != (len(seq) - 1) and seq[pos + 1] != '-' and seq[pos + 1] != 'p' and seq[pos + 1] != 'P' and seq[pos + 1] != 'z' and seq[pos + 1] != 'Z':
        new_seq = seq[:pos] + '-' + seq[pos:]
        seq = new_seq
    return seq

def demethylate(seq):
    if '-' in seq:
        pos = find_methylation_in_seq(seq)
        if not pos:
            return seq
        new_seq = seq[:pos] + seq[pos + 1:]
        seq = new_seq
    return seq


def smiles_from_seq(seq, cyclize):
    """Calculates the smiles of a given peptide dendrimer sequence

    Arguments:
        seq {string} -- peptide dendrimer sequence
    Returns:
        string -- molecule_smile - SMILES of the peptide
    """

    #seq = seq.replace("-z","z").replace("-Z","Z").replace("-p","p").replace("-P","P")

    gs, bs, terminal, capping = split_seq_components(seq)

    # modifies the Cterminal
    if terminal:
        molecule = rdmolfiles.MolFromSmiles(T_SMILES[terminal[0]])
    else:
        molecule = ''


    if cyclize and bs:
        print('dendrimer, cyclization not possible, branching unit will not be considered')

    if cyclize:
        for gen in gs:
            metbond = False
            for aa in gen:
                if aa == 'X':
                    continue
                if aa == '-':
                    metbond = True
                    continue
                if molecule == '':
                    molecule = rdmolfiles.MolFromSmiles(AA_SMILES[aa])
                else:
                    molecule = utils.connect_mol(molecule, rdmolfiles.MolFromSmiles(AA_SMILES[aa]), metbond)
                    if metbond:
                        metbond = False
    else:
        # creates the dendrimer structure
        for gen in gs:
            metbond = False
            for aa in gen:
                if aa == '-':
                    metbond = True
                    continue
                if molecule == '':
                    molecule = rdmolfiles.MolFromSmiles(AA_SMILES[aa])
                else:
                    molecule = utils.connect_mol(molecule, rdmolfiles.MolFromSmiles(AA_SMILES[aa]), metbond)
                    if metbond:
                        metbond = False

            if bs:
                if bs[0] == '-':
                    metbond = True
                    bs.pop(0)
                if molecule == '':
                    molecule = rdmolfiles.MolFromSmiles(B_SMILES[bs[0]])
                else:
                    molecule = utils.connect_mol(molecule, rdmolfiles.MolFromSmiles(B_SMILES[bs[0]]), metbond)
                    if metbond:
                        metbond = False
                bs.pop(0)

    # adds capping to the N-terminal (the called clip function is different, cause the listed smiles 
    # for the capping are already without OH, it is not necessary removing any atom after foming the new bond)
    
    if molecule == '':
        smiles = ''
        return smiles, seq
    
    if capping:
        molecule = utils.attach_capping(molecule, rdmolfiles.MolFromSmiles(C_SMILES[capping[0]]))

    if cyclize:
        if is_cyclic(seq):
            cy = 1
        else:
            cy = 0
        molecule = utils.cyclize(molecule, cy)

    # clean the smile from all the tags
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)

    molecule_smile = rdmolfiles.MolToSmiles(molecule, isomericSmiles=True).replace('[N]', 'N').replace('[C]', 'C')
    return molecule_smile, seq


def remove_SS_cyclization(seq):
    if 'Ä' in seq:
        seq = seq.replace('Ä', '')
    if 'Ö' in seq:
        seq = seq.replace('Ö', '')
    if 'Ü' in seq:
        seq = seq.replace('Ü', '')
    return seq    
    
def remove_methylation_before_proline(seq):
    # not methylation before proline
    return seq.replace("-Z", "Z").replace("-z", "z").replace("-P", "P").replace("-p", "p")

def remove_methylation_at_end_and_begining(seq):
    # no methylation when sequence ends
    if not is_empty(seq):
        if seq[-1] == "-":
            seq = seq[:-1]
    if not is_empty(seq):
        if seq[0] == "-":
            seq = seq[1:]
    if not is_empty(seq):
        if seq[0] == "X" and seq[1] == "-":
            seq = seq[0] + seq[2:]
    return seq

def swapcy(seq):
    """insertion of two ativated cys at head to tail position

    Arguments:
        seq {string} -- peptide seq

    Returns:
        string -- S-S cyclized peptide seq
    """
    if is_cyclic(seq):
        act_cys = 'Ä'
        if 'Ä' in seq:
            act_cys = 'Ö'
            if 'Ö' in seq:
                act_cys = 'Ü'
                if 'Ü' in seq:
                    return seq

        new_seq = act_cys + seq[1:] + act_cys

        return new_seq
    else:
        return seq


def break_SS(seq):
    """inactivation of all cys

    Arguments:
        seq {string} -- peptide seq

    Returns:
        string -- S-S cyclized peptide seq
    """

    act_cys = 'Ü'
    if 'Ü' not in seq:
        act_cys = 'Ö'
        if 'Ö' not in seq:
            act_cys = 'Ä'
            if 'Ä' not in seq:
                return seq

    # seq.replace('Ä', 'C')
    # seq.replace('Ü', 'C')
    # seq.replace('Ö' ,'C')
    seq.replace(act_cys, '')

    return seq
