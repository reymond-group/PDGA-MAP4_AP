import random

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdmolops

from . import sequence

def cyclize(mol, cy):
    """it is connecting cyclizing the given molecule

    Arguments:
        mol {rdKit mol object} -- molecule to be cyclized
        cy {int} -- 1=yes, 0=no cyclazation

    Returns:
        mols {list of rdKit mol objects} -- possible cyclazation
    """
    count = 0

    # detects all the N terminals in mol
    for atom in mol.GetAtoms():
        if atom.GetSmarts() == '[N:2]' or atom.GetSmarts() == '[NH2:2]' or atom.GetSmarts() == '[NH:2]':
            count += 1
            atom.SetProp('Nterm', 'True')
        else:
            atom.SetProp('Nterm', 'False')

    # detects all the C terminals in mol (it should be one)
    for atom in mol.GetAtoms():
        if atom.GetSmarts() == '[C:1]' or atom.GetSmarts() == '[CH:1]':
            atom.SetProp('Cterm', 'True')
        else:
            atom.SetProp('Cterm', 'False')

    # detects all the S terminals in mol

    for atom in mol.GetAtoms():
        if atom.GetSmarts() == '[S:1]':
            atom.SetProp('Sact1', 'True')
        else:
            atom.SetProp('Sact1', 'False')

    for atom in mol.GetAtoms():
        if atom.GetSmarts() == '[S:2]':
            atom.SetProp('Sact2', 'True')
        else:
            atom.SetProp('Sact2', 'False')

    for atom in mol.GetAtoms():
        if atom.GetSmarts() == '[S:3]':
            atom.SetProp('Sact3', 'True')
        else:
            atom.SetProp('Sact3', 'False')

    Nterm = []
    Cterm = []
    Sact1 = []
    Sact2 = []
    Sact3 = []

    # saves active Cysteins postions:
    for atom in mol.GetAtoms():
        if atom.GetProp('Sact1') == 'True':
            Sact1.append(atom.GetIdx())

    # saves active Cysteins 2 postions:
    for atom in mol.GetAtoms():
        if atom.GetProp('Sact2') == 'True':
            Sact2.append(atom.GetIdx())

    # saves active Cysteins 3 postions:
    for atom in mol.GetAtoms():
        if atom.GetProp('Sact3') == 'True':
            Sact3.append(atom.GetIdx())

    # creates the S-S bond (in the current version only two 'active' Cys, this codo picks two random anyway):
    while len(Sact1) >= 2:
        edmol = rdchem.EditableMol(mol)
        pos = list(range(len(Sact1)))
        x = np.random.choice(pos, 1)[0]
        pos.remove(x)
        y = np.random.choice(pos, 1)[0]
        a = Sact1[x]
        b = Sact1[y]
        edmol.AddBond(a, b, order=Chem.rdchem.BondType.SINGLE)
        mol = edmol.GetMol()
        mol.GetAtomWithIdx(a).SetProp('Sact1', 'False')
        mol.GetAtomWithIdx(b).SetProp('Sact1', 'False')
        mol.GetAtomWithIdx(a).SetAtomMapNum(0)
        mol.GetAtomWithIdx(b).SetAtomMapNum(0)
        Sact1.remove(a)
        Sact1.remove(b)

    while len(Sact2) >= 2:
        edmol = rdchem.EditableMol(mol)
        pos = list(range(len(Sact2)))
        x = np.random.choice(pos, 1)[0]
        pos.remove(x)
        y = np.random.choice(pos, 1)[0]
        a = Sact2[x]
        b = Sact2[y]
        edmol.AddBond(a, b, order=Chem.rdchem.BondType.SINGLE)
        mol = edmol.GetMol()
        mol.GetAtomWithIdx(a).SetProp('Sact2', 'False')
        mol.GetAtomWithIdx(b).SetProp('Sact2', 'False')
        mol.GetAtomWithIdx(a).SetAtomMapNum(0)
        mol.GetAtomWithIdx(b).SetAtomMapNum(0)
        Sact2.remove(a)
        Sact2.remove(b)

    while len(Sact3) >= 2:
        edmol = rdchem.EditableMol(mol)
        pos = list(range(len(Sact3)))
        x = np.random.choice(pos, 1)[0]
        pos.remove(x)
        y = np.random.choice(pos, 1)[0]
        a = Sact3[x]
        b = Sact3[y]
        edmol.AddBond(a, b, order=Chem.rdchem.BondType.SINGLE)
        mol = edmol.GetMol()
        mol.GetAtomWithIdx(a).SetProp('Sact3', 'False')
        mol.GetAtomWithIdx(b).SetProp('Sact3', 'False')
        mol.GetAtomWithIdx(a).SetAtomMapNum(0)
        mol.GetAtomWithIdx(b).SetAtomMapNum(0)
        Sact3.remove(a)
        Sact3.remove(b)

    # saves active C and N terminals postions:
    for atom in mol.GetAtoms():
        if atom.GetProp('Nterm') == 'True':
            Nterm.append(atom.GetIdx())
        if atom.GetProp('Cterm') == 'True':
            Cterm.append(atom.GetIdx())

    if cy == 1:
        edmol = rdchem.EditableMol(mol)

        # creates the amide bond
        edmol.AddBond(Nterm[0], Cterm[0], order=Chem.rdchem.BondType.SINGLE)
        edmol.RemoveAtom(Cterm[0] + 1)

        mol = edmol.GetMol()

        # removes tags and lables form the atoms which reacted
        mol.GetAtomWithIdx(Nterm[0]).SetProp('Nterm', 'False')
        mol.GetAtomWithIdx(Cterm[0]).SetProp('Cterm', 'False')
        mol.GetAtomWithIdx(Nterm[0]).SetAtomMapNum(0)
        mol.GetAtomWithIdx(Cterm[0]).SetAtomMapNum(0)

    return mol


def attach_capping(mol1, mol2):
    """it is connecting all Nterminals with the desired capping

    Arguments:
        mol1 {rdKit mol object} -- first molecule to be connected
        mol2 {rdKit mol object} -- second molecule to be connected - chosen N-capping

    Returns:
        rdKit mol object -- mol1 updated (connected with mol2, one or more)
    """

    count = 0

    # detects all the N terminals in mol1
    for atom in mol1.GetAtoms():
        atom.SetProp('Cterm', 'False')
        if atom.GetSmarts() == '[N:2]' or atom.GetSmarts() == '[NH2:2]' or atom.GetSmarts() == '[NH:2]':
            count += 1
            atom.SetProp('Nterm', 'True')
        else:
            atom.SetProp('Nterm', 'False')

    # detects all the C terminals in mol2 (it should be one)
    for atom in mol2.GetAtoms():
        atom.SetProp('Nterm', 'False')
        if atom.GetSmarts() == '[C:1]' or atom.GetSmarts() == '[CH:1]':
            atom.SetProp('Cterm', 'True')
        else:
            atom.SetProp('Cterm', 'False')

    # mol2 is addes to all the N terminal of mol1
    for i in range(count):
        combo = rdmolops.CombineMols(mol1, mol2)
        Nterm = []
        Cterm = []

        # saves in two different lists the index of the atoms which has to be connected
        for atom in combo.GetAtoms():
            if atom.GetProp('Nterm') == 'True':
                Nterm.append(atom.GetIdx())
            if atom.GetProp('Cterm') == 'True':
                Cterm.append(atom.GetIdx())

        # creates the amide bond
        edcombo = rdchem.EditableMol(combo)
        edcombo.AddBond(Nterm[0], Cterm[0], order=Chem.rdchem.BondType.SINGLE)
        clippedMol = edcombo.GetMol()

        # removes tags and lables form the atoms which reacted
        clippedMol.GetAtomWithIdx(Nterm[0]).SetProp('Nterm', 'False')
        clippedMol.GetAtomWithIdx(Cterm[0]).SetProp('Cterm', 'False')
        clippedMol.GetAtomWithIdx(Nterm[0]).SetAtomMapNum(0)
        clippedMol.GetAtomWithIdx(Cterm[0]).SetAtomMapNum(0)
        # uptades the 'core' molecule
        mol1 = clippedMol

    return mol1



def remove_duplicates(gen):
    """Removes duplicates

    Arguments:
        gen {list} -- sequences

    Returns:
        list -- unique list of sequences
    """

    gen_u = []
    for seq in gen:
        if seq not in gen_u:
            gen_u.append(seq)
    return gen_u


def mating(parents):
    """splits the parents in half and join them giving a child

    Arguments:
        parents {list of strings} -- parents

    Returns:
        string -- child
    """

    parent1 = parents[0]
    parent2 = parents[1]
    half1 = parent1[:random.randint(int(round(len(parent1) / 2, 0)) - 1, int(round(len(parent1) / 2, 0)) + 1)]
    half2 = parent2[random.randint(int(round(len(parent2) / 2, 0)) - 1, int(round(len(parent2) / 2, 0)) + 1):]
    child = half1 + half2

    child = sequence.sanitize_sequence(child)

    return child



def set_seed(seed):
    """set seed for random

    Arguments:
        seed {int} -- sed for random
    """

    random.seed(int(seed))
    np.random.seed(int(seed))


def print_time(time):
    """print running time
    """

    hours, rem = divmod(time, 3600)
    minutes, seconds = divmod(rem, 60)
    print('Time {:0>2}:{:0>2}:{:0>2}'.format(int(hours), int(minutes), int(seconds)))


def write_progress(path, dist_dict, gen_n, jd_av, jd_min ):
    """add gen number, gen sequences and its jd av and min
    
    """

    gen_temp = []
    gen = list(dist_dict.keys())
    for seq in gen:
        gen_temp.append(sequence.reinterprete(seq))
    gen = ';'.join(map(str, gen_temp))
    with open(path , 'a') as outFile:
        outFile.write(str(gen_n) + ' ' + gen + ' ' + str(jd_av) + ' ' + str(jd_min) + '\n')

def write_results(path, smiles, seq, map4, jd):
    """if jd from query is smaller than similarity treshold
        (class variable), adds seq to results
    
    """
    with open(path, 'a') as outFile:
        outFile.write(smiles + ' ' + sequence.reinterprete(seq) + ' ' + str(round(jd,3)) + '\n')


def connect_mol(mol1, mol2, metbond):
        """it is connecting all Nterminals of mol1 with the Cterminal 
        of the maximum possible number of mol2s
    
        Arguments:
            mol1 {rdKit mol object} -- first molecule to be connected
            mol2 {rdKit mol object} -- second molecule to be connected

        Returns:
            rdKit mol object -- mol1 updated (connected with mol2, one or more)
        """
        count = 0

        # detects all the N terminals in mol1
        for atom in mol1.GetAtoms():
            atom.SetProp('Cterm', 'False')
            atom.SetProp('methyl', 'False')
            if atom.GetSmarts() == '[N:2]' or atom.GetSmarts() == '[NH2:2]' or atom.GetSmarts() == '[NH:2]':
                count += 1
                atom.SetProp('Nterm', 'True')
            else:
                atom.SetProp('Nterm', 'False')

        # detects all the C terminals in mol2 (it should be one)
        for atom in mol2.GetAtoms():
            atom.SetProp('Nterm', 'False')
            atom.SetProp('methyl', 'False')
            if atom.GetSmarts() == '[C:1]' or atom.GetSmarts() == '[CH:1]':
                atom.SetProp('Cterm', 'True')
            else:
                atom.SetProp('Cterm', 'False')

        # mol2 is addes to all the N terminal of mol1
        for i in range(count):
            combo = rdmolops.CombineMols(mol1, mol2)
            Nterm = []
            Cterm = []

            # saves in two different lists the index of the atoms which has to be connected
            for atom in combo.GetAtoms():
                if atom.GetProp('Nterm') == 'True':
                    Nterm.append(atom.GetIdx())
                if atom.GetProp('Cterm') == 'True':
                    Cterm.append(atom.GetIdx())

            # creates the amide bond
            edcombo = rdchem.EditableMol(combo)
            edcombo.AddBond(Nterm[0], Cterm[0], order=Chem.rdchem.BondType.SINGLE)
            edcombo.RemoveAtom(Cterm[0] + 1)
            clippedMol = edcombo.GetMol()

            # removes tags and lables form c term atoms which reacted
            clippedMol.GetAtomWithIdx(Cterm[0]).SetProp('Cterm', 'False')
            clippedMol.GetAtomWithIdx(Cterm[0]).SetAtomMapNum(0)

            # methylates amide bond
            if metbond == True:
                Nterm = []
                Met = []
                methyl = rdmolfiles.MolFromSmiles('[C:4]')
                for atom in methyl.GetAtoms():
                    atom.SetProp('methyl', 'True')
                    atom.SetProp('Nterm', 'False')
                    atom.SetProp('Cterm', 'False')
                metcombo = rdmolops.CombineMols(clippedMol, methyl)
                for atom in metcombo.GetAtoms():
                    if atom.GetProp('Nterm') == 'True':
                        Nterm.append(atom.GetIdx())
                    if atom.GetProp('methyl') == 'True':
                        Met.append(atom.GetIdx())
                metedcombo = rdchem.EditableMol(metcombo)
                metedcombo.AddBond(Nterm[0], Met[0], order=Chem.rdchem.BondType.SINGLE)
                clippedMol = metedcombo.GetMol()
                clippedMol.GetAtomWithIdx(Met[0]).SetProp('methyl', 'False')
                clippedMol.GetAtomWithIdx(Met[0]).SetAtomMapNum(0)

            # removes tags and lables form the atoms which reacted
            clippedMol.GetAtomWithIdx(Nterm[0]).SetProp('Nterm', 'False')
            clippedMol.GetAtomWithIdx(Nterm[0]).SetAtomMapNum(0)

            # uptades the 'core' molecule
            mol1 = clippedMol
        
        return mol1


def random_subset(stuff, n, seed=None):
    rng = np.random.default_rng(seed=seed)
    return rng.choice(stuff, n, replace=False)