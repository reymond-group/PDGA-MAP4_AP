import random 
from . import sequence

class SequenceGenerator:
    def __init__(self, verbose=False):
            # variables for random generation of dendrimers
            self.AA4rndm = ['O', 'Z', 'R', 'H', 'K', 'E', 'S', 'T', 'N', 'Q', 'G', 'P', 'A', 'V', 'I', 'L', 'F', 'Y', 'W', 'C', 'D',
                    'M',
                    '!', '?', '=', '%', '$', '@', '#', '', '', '', '', '', '', '',
                    '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
            self.B4rndm = ['1', '2', '3', '4', '']
            self.CTrndm = ['+', '']
            self.NTrndm = ['&', '']

            self.max_aa_no = 5
            self.max_gen_no = 3

            self.verbose = verbose
    
    def __repr__(self):
        return "SequenceGenerator(AA4rndm={0},B4rndm={1},CTrndm={2},NTrndm={3},max_aa_no={4},max_gen_no={5},verbose={6})".format(self.AA4rndm, self.B4rndm, self.CTrndm, self.NTrndm, self.max_aa_no, self.max_gen_no, self.verbose)


    def generate(self):
        """Generates random implicit sequences of max "max_gen_no" generation dendrimers
            with max "max_aa_no" AA in each generation, picking from AA4random, B4random
            (probability of position to be empty intrinsic in these lists). 
        
        Returns:
            string -- implicit sequence of a random dendrimer
        """

        new_random_seq = random.choice(self.CTrndm)
        aa_count = 0

        while aa_count < self.max_aa_no:
            new_random_seq += random.choice(self.AA4rndm)
            aa_count += 1
        gen_count = 0
        while gen_count < self.max_gen_no:
            if new_random_seq != '':
                new_random_seq += random.choice(self.B4rndm)
            aa_count = 0
            while aa_count < self.max_aa_no:
                new_random_seq += random.choice(self.AA4rndm)
                aa_count += 1
            gen_count += 1
        new_random_seq += random.choice(self.NTrndm)

        return new_random_seq[::-1]

    def exclude_buildingblocks(self, bb_to_ex):
        """Excludes the given building blocks
        
        Arguments:
            bb_to_ex {list} -- building blocks to exclude
        """

        for bb in bb_to_ex:
            if bb in sequence.interprete_dict.keys():
                element = sequence.interprete(bb)
                if element in sequence.AA:
                    self.exclude_aminoacids(element)
                elif element in sequence.B:
                    self.exclude_branching(element)
                elif element in sequence.CT:
                    self.exclude_C_terminal(element)
                elif element in sequence.NT:
                    self.exclude_N_capping(element)                
                else:
                    print("can't exclude ", bb)
            else:
                print("can't exclude ", bb)


    def exclude_aminoacids(self, aa_to_ex):
        """Excludes the given aminoacids
        
        Arguments:
            aa_to_ex {list} -- aminoacids to exclude
        """

        for element in aa_to_ex:
            sequence.AA.remove(element)
            self.AA4rndm.remove(element)
            self.AA4rndm.remove('')

        if not sequence.AA:
            sequence.AA.append('')
        if not self.AA4rndm:
            self.AA4rndm.append('')

        if self.verbose:
            print('The GA is using aminoacids:', sequence.AA)

    def exclude_branching(self, bb_to_ex):
        """Excludes the given branching units
        
        Arguments:
            bb_to_ex {list} -- branching units to exclude
        """

        if self.porpouse != 'cyclic':

            for element in bb_to_ex:
                self.B.remove(element)
                self.B4rndm.remove(element)

            if not self.B:
                self.B.append('')
            if not self.B4rndm:
                self.B4rndm.append('')

            if self.verbose:
                print('The GA is using branching units:', self.B)

    def exclude_C_terminal(self, t_to_ex):
        """Excludes the given C terminal modifications
        
        Arguments:
            t_to_ex {list} -- C terminal modifications to exclude
        """

        for element in t_to_ex:
            self.CT.remove(element)
            self.CTrndm.remove(element)
            self.CTrndm.remove('')

        if not self.CT:
            self.CT.append('')
        if not self.CTrndm:
            self.CTrndm.append('')

        if self.verbose:
            print('The GA is using C terminal mod:', self.CT)

    def exclude_N_capping(self, c_to_ex):
        """Excludes the given N terminal capping
        
        Arguments:
            c_to_ex {list} -- N terminal capping to exclude
        """

        for element in c_to_ex:
            self.NT.remove(element)
            self.NTrndm.remove(element)
            self.NTrndm.remove('')

        if not self.NT:
            self.NT.append('')
        if not self.NTrndm:
            self.NTrndm.append('')

        if self.verbose:
            print('The GA is using N capping mod:', self.NT)
