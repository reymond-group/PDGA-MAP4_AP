from typing import Sequence
from . import sequence
import numpy as np 


def delete(seq, mut_n=1):
    """Performs n (mut_n, class variable) deletion 

    Arguments:
        seq {string} -- seq to be mutated

    Returns:
        list -- mutations
    """

    mutations = []
    for i in range(mut_n):
        pos = sequence.find_anything_in_seq(seq)
        if not pos:
            continue
        seq1 = seq[:pos]
        seq2 = seq[pos + 1:]
        new_seq = seq1 + seq2
        mutations.append(new_seq)

    return mutations


def insert(seq, type_insert, mut_n=1):
        """Performs n (mut_n, class variable) random point insertions. 
        If type insert is 'aa' the new element will be an aminoacid.
        If type insert is 'b' the new element will be a branching unit.
    
        Arguments:
            seq {string} -- seq to be mutated

        Returns:
            list -- mutations
        """

        mutations = []
        for i in range(mut_n):
            pos = sequence.find_anything_in_seq(seq)
            if not pos:
                continue
            if type_insert == 'aa':
                new_element = np.random.choice(sequence.AA, 1)
            elif type_insert == 'b':
                new_element = np.random.choice(sequence.B, 1)
            else:
                raise ValueError("not valid type, type has to be \"aa\" or \"b\"")

            seq1 = seq[:pos]
            seq2 = seq[pos:]
            seq = seq1 + new_element[0] + seq2
            mutations.append(seq)

        return mutations

def mutate_aa(seq, mut_n=1):
    """Performs n (mut_n, class variable) random point mutation

    Arguments:
        seq {string} -- seq to be mutated

    Returns:
        list -- mutations
    """
    mutations = []

    for i in range(mut_n):
        aa_pos = sequence.find_aminoacid_in_seq(seq)
        if not aa_pos:
            continue
        seq1 = seq[:aa_pos]
        seq2 = seq[aa_pos + 1:]
        aa_new = np.random.choice(sequence.AA, 1)
        seq = seq1 + aa_new[0] + seq2
        mutations.append(seq)

    return mutations

def move_branching_point(seq, pos=+1, mut_n=1, verbose=False):
    """Performs n (mut_n, class variable) random point mutation

    Arguments:
        seq {string} -- seq to be mutated
        pos {integer} -- position to move the branching unit, positive for right, negative for left

    Returns:
        list -- mutations
    """

    mutations = []
    for i in range(mut_n):
        b_pos = sequence.find_branching_in_seq(seq)
        if not b_pos:
            continue
        b = seq[b_pos]
        if 0 <= b_pos + pos < len(seq):
            if seq[b_pos + pos] in sequence.CT or seq[b_pos + pos] in sequence.NT:
                mutations.append(seq)
                if verbose:
                    print(seq + ' Terminal found, could not move ' + b + ' {}'.format(pos))
                continue
            else:
                seqd = seq[:b_pos] + seq[b_pos + 1:]
                seq1 = seqd[:b_pos + pos]
                seq2 = seqd[b_pos + pos:]
                seq = seq1 + b + seq2
                mutations.append(seq)
        else:
            mutations.append(seq)

    return mutations


def mutate_branching_point(seq, mut_n=1):
    """Performs n (mut_n, class variable) random point mutation

    Arguments:
        seq {string} -- seq to be mutated

    Returns:
        list -- mutations
    """
    mutations = []

    for i in range(mut_n):
        b_pos = sequence.find_branching_in_seq(seq)
        if not b_pos:
            continue
        seq1 = seq[:b_pos]
        seq2 = seq[b_pos + 1:]
        b_new = np.random.choice(sequence.B, 1)
        seq = seq1 + b_new[0] + seq2
        mutations.append(seq)

    return mutations



def mutation_deletion(gen, mut_rate):
    gen_tmp = []
    seq_deletion = np.random.choice(gen, int(round(len(gen) * mut_rate, 0)), replace=False)
    for seq in gen:
        if seq in seq_deletion:
            aa_pos, b_pos, met_pos, all_pos = sequence.find_aa_b_pos(seq)
            if all_pos:
                seq_delete = delete(seq)
                for s in seq_delete:
                    gen_tmp.append(s)
            else:
                gen_tmp.append(seq)
        else:
            gen_tmp.append(seq)
    return gen_tmp


def mutation_insertion_aa(gen, mut_rate):
    gen_tmp = []
    seq_insertion_aa = np.random.choice(gen, int(round(len(gen) * mut_rate, 0)), replace=False)
    for seq in gen: 
        if seq in seq_insertion_aa:
            aa_pos, b_pos, met_pos, all_pos = sequence.find_aa_b_pos(seq)
            if all_pos:
                seq_insert = insert(seq, 'aa')
                for s in seq_insert:
                    gen_tmp.append(s)
            else:
                gen_tmp.append(seq)
        else:
            gen_tmp.append(seq)
    return gen_tmp



def mutation_insertion_branch(gen, b_insert_rate):
    gen_tmp = []
    # to avoid incontrolled progressivly growth of the sequences,
    # mantain b_insert_rate (class variable, default = 0.1) low
    seq_insertion_b = np.random.choice(gen, int(round(len(gen) * b_insert_rate, 0)), replace=False)
    for seq in gen:
        if seq in seq_insertion_b:
            aa_pos, b_pos, met_pos, all_pos = sequence.find_aa_b_pos(seq)
            if all_pos:
                seq_insert = insert(seq, 'b')
                for s in seq_insert:
                    gen_tmp.append(s)                    
            else:
                gen_tmp.append(seq)
        else:
            gen_tmp.append(seq)
    return gen_tmp


def mutation_mutate_aa(gen, mut_rate):
    gen_tmp = []      
    seqs_mutate_aa = np.random.choice(gen, int(round(len(gen) * mut_rate, 0)), replace=False)
    for seq in gen:
        if seq in seqs_mutate_aa:
            aa_pos, b_pos, met_pos, all_pos = sequence.find_aa_b_pos(seq)
            if aa_pos:
                seq_mutate = mutate_aa(seq)
                for s in seq_mutate:
                    gen_tmp.append(s)
            else:
                gen_tmp.append(seq)
        else:
            gen_tmp.append(seq)
    return gen_tmp

def mutate_move_branch(gen, mut_rate, position):
    gen_tmp = []      
    seq_move_b = np.random.choice(gen, int(round(len(gen) * mut_rate, 0)), replace=False)
    for seq in gen:
        if seq in seq_move_b:
            aa_pos, b_pos, met_pos, all_pos = sequence.find_aa_b_pos(seq)
            if b_pos:
                seq_mutate = move_branching_point(seq, position)
                for s in seq_mutate:
                    gen_tmp.append(s)
            else:
                gen_tmp.append(seq)
        else:
            gen_tmp.append(seq)
    return gen_tmp

def mutate_move_branch_right(gen, mut_rate):
    return mutate_move_branch(gen,mut_rate, +1)

def mutate_move_branch_left(gen, mut_rate):
    return mutate_move_branch(gen,mut_rate, -1)


def mutate_b(gen, mut_rate):
    gen_tmp = []      
    seqs_mutate_b = np.random.choice(gen, int(round(len(gen) * mut_rate, 0)), replace=False)
    for seq in gen:
        if seq in seqs_mutate_b:
            aa_pos, b_pos, met_pos, all_pos = sequence.find_aa_b_pos(seq)
            if b_pos:
                seq_mutate = mutate_branching_point(seq)
                for s in seq_mutate:
                    gen_tmp.append(s)
            else:
                gen_tmp.append(seq)           
        else:
            gen_tmp.append(seq)
    return gen_tmp


def mutate_c(gen, mut_rate):
    gen_tmp = []      
    seqs_mutate_c = np.random.choice(gen, int(round(len(gen) * mut_rate, 0)), replace=False)
    for seq in gen:
        if seq in seqs_mutate_c:
            
            if sequence.is_cyclic(seq):
                gen_tmp.append(seq)
                continue

            c = np.random.choice(sequence.NT, 1, replace=False)[0]
            if len(seq) > 2 and seq[0] in sequence.NT:
                seq_tmp = seq[1:]
                new_seq = c + seq_tmp
                gen_tmp.append(new_seq)
            else:
                gen_tmp.append(seq)
        else:
            gen_tmp.append(seq)
    return gen_tmp


def mutate_t(gen, mut_rate):
    gen_tmp = []      
    seqs_mutate_t = np.random.choice(gen, int(round(len(gen) * mut_rate, 0)), replace=False)
    for seq in gen:
        if seq in seqs_mutate_t:

            if sequence.is_cyclic(seq):
                gen_tmp.append(seq)
                continue

            t = np.random.choice(sequence.CT, 1, replace=False)[0]
            if len(seq) > 2 and seq[-1] in sequence.CT:
                seq_tmp = seq[:-1]
                new_seq = seq_tmp + t
                gen_tmp.append(new_seq)
            else:
                gen_tmp.append(seq)
        else:
            gen_tmp.append(seq)
    return gen_tmp


def mutate_methylation(gen, mut_rate):
    gen_tmp = []      
    seqs_methylate = np.random.choice(gen, int(round(len(gen) * mut_rate, 0)), replace=False)
    for seq in gen:
        if seq in seqs_methylate:
            gen_tmp.append(sequence.methylate(seq))
        else:
            gen_tmp.append(seq)
    return gen_tmp

def mutate_demethylation(gen, mut_rate):
    gen_tmp = []      
    seqs_demethylate = np.random.choice(gen, int(round(len(gen) * mut_rate, 0)), replace=False)
    for seq in gen:
        if seq in seqs_demethylate:
            gen_tmp.append(sequence.demethylate(seq))
        else:
            gen_tmp.append(seq)        
    return gen_tmp


def mutate_break_ss(gen, mut_rate):
    gen_tmp = []      
    seqs_inact_cys = np.random.choice(gen, int(round(len(gen) * mut_rate, 0)), replace=False)
    for seq in gen:
        if seq in seqs_inact_cys:
            gen_tmp.append(sequence.break_SS(seq))
        else:
            gen_tmp.append(seq)
    return gen_tmp

def mutate_make_ss(gen, mut_rate):
    gen_tmp = []      
    seqs_act_cys = np.random.choice(gen, int(round(len(gen) * mut_rate, 0)), replace=False)
    for seq in gen:
        if seq in seqs_act_cys:
            act_seq = sequence.form_SS(seq)
            gen_tmp.append(act_seq)
        else:
            gen_tmp.append(seq)
    return gen_tmp

def mutate_linearize(gen, mut_rate):
    gen_tmp = []      
    seqs_lin = np.random.choice(gen, int(round(len(gen) * mut_rate, 0)), replace=False)
    for seq in gen:
        if seq in seqs_lin:
            new_seq = sequence.linearize(seq)
            gen_tmp.append(new_seq)
        else:
            gen_tmp.append(seq)
    return gen_tmp

def mutate_cyclize(gen, mut_rate):
    gen_tmp = []      
    seqs_cy = np.random.choice(gen, int(round(len(gen) * mut_rate, 0)), replace=False)
    for seq in gen:
        if seq in seqs_cy:
            new_seq = sequence.add_cycle(seq)
            gen_tmp.append(new_seq)
        else:
            gen_tmp.append(seq)
    return  gen_tmp


def mutate_head_tail_ss(gen, mut_rate):
    gen_tmp = []      
    # swap head-to-tail with S-S
    seqs_swap = np.random.choice(gen, int(round(len(gen) * mut_rate, 0)), replace=False)
    for seq in gen:
        if seq in seqs_swap:
            new_seq = sequence.swapcy(seq)
            gen_tmp.append(new_seq)
        else:
            gen_tmp.append(seq)
    return gen_tmp




linear_mutations = { 
    "deletion": mutation_deletion,
    "insertion_aa" : mutation_insertion_aa,
    "insertion_branch" : mutation_insertion_branch,
    "mutation_aa": mutation_mutate_aa,
    "branch_right": mutate_move_branch_right,
    "branch_left": mutate_move_branch_left,
    "mutate_branch": mutate_b,
    "mutate_cterm": mutate_c,
    "mutate_nterm": mutate_t,
    "add_methylation" : mutate_methylation,
    "remove_methylation" : mutate_demethylation
    }

cyclic_mutations = { 
    "deletion": mutation_deletion,
    "insertion_aa" : mutation_insertion_aa,
    "mutation_aa": mutation_mutate_aa,
    "mutate_cterm": mutate_c,
    "mutate_nterm": mutate_t,
    "break_ss" : mutate_break_ss,
    "make_ss": mutate_make_ss,
    "linearize": mutate_linearize,
    "cyclize": mutate_cyclize,
    "head_tail_ss": mutate_head_tail_ss,
    "add_methylation" : mutate_methylation,
    "remove_methylation" : mutate_demethylation
    }


def mutate(gen, cyclic, mut_rate=1, b_insert_rate=0.1, methyl=False):
    """Mutates the given generation 

    Arguments:
        gen {list} -- sequences


    Returns:
        [list] -- mutated generation
    """

    if cyclic:
        mutations = cyclic_mutations
    else:
        mutations = linear_mutations
    gen_tmp = []

    for mutation in np.random.choice(list(mutations.keys()), 1, replace=False):
        if mutation == "add_methylation" or mutation == "remove_methylation":
            if methyl:
                gen_tmp.extend(mutations[mutation](gen, mut_rate))
        else:
            if mutation == "insertion_branch":
                gen_tmp.extend(mutations[mutation](gen, b_insert_rate))
            else:
                gen_tmp.extend(mutations[mutation](gen, mut_rate))


    if gen_tmp == []:
        gen_tmp = gen
    gen_new = []

    for seq in gen_tmp:
        seq = sequence.sanitize_sequence(seq)
        if not sequence.is_empty(seq):
            gen_new.append(seq)

    return gen_new