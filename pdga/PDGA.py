import os
import time

import numpy as np
from rdkit import Chem

from . import mutations
from . import sequence
from . import utils
from .sequence_random_generator import SequenceGenerator


class PDGA:
    def __init__(self, pop_size, mut_rate, gen_gap, query, sim_treshold, porpouse, folder,\
         fingerprintfn, distancefn, query_name, is_peptide_sequence=True, methyl=False, verbose=False, seed=None):
        
        self.pop_size = int(pop_size)
        self.mut_rate = float(mut_rate)
        self.gen_gap = float(gen_gap)
        self.porpouse = porpouse


        if not os.path.exists(folder):
            os.makedirs(folder)
        self.folder = folder
        if is_peptide_sequence:
            self.query = sequence.interprete(query)
        else:
            self.query = query

        
        self.fingerprintfn = fingerprintfn
        self.distancefn = distancefn

        proc_seqs, fps, smiles = self.calc_fingerprints([self.query], is_peptide_sequence=is_peptide_sequence)
        self.query_fp = fps[0]

        self.sim_treshold = sim_treshold

        self.output_path = '{}/{}_generations'.format(folder, query_name)
        self.results_path = '{}/{}_results'.format(folder, query_name)

        self.timelimit_seconds = None

        self.mut_n = 1
        self.b_insert_rate = 0.1
        self.selec_strategy = 'Elitist'
        self.rndm_newgen_fract = 10

        self.gen_n = 0

        self.sequence_rng = SequenceGenerator()
        if self.porpouse == 'linear' or self.porpouse == 'cyclic':
            sequence.B = ['']
            self.sequence_rng.B4rndm = ['']

        self._methyl = methyl
        self.verbose = verbose


        if seed is None:
            seed = np.random.randint(100000000)
            if verbose:
                print("Use seed:",seed)
        self.seed = seed
        utils.set_seed(seed=seed)


    def rndm_gen(self):
        """Creates a generation of "pop_size" random dendrimers        
        Returns:
           list -- generation of "pop_size" random dendrimers
        """

        gen = []
        while len(gen) < self.pop_size:
            gen.append(self.sequence_rng.generate())
        return gen


    def calc_fingerprints(self, seqs, is_peptide_sequence=True):
        """Calculates the map4 for the given values
        
        Arguments:
            seqs {list} -- peptide sequence list
        
        Returns:
            precessed sesq, fps, smiles, props {lists} -- processed sequences and the relative map4
        """

        proc_seqs = []
        smiles = []
        mol_list = []

        
        for seq in seqs:
            if seq == '':
                continue

            if is_peptide_sequence:
                smi, seq = sequence.smiles_from_seq(seq, cyclize = self.porpouse == 'cyclic')
            else:
                smi = seq

            if smi == '':
                continue

            mol = Chem.MolFromSmiles(smi) 
            if mol:
                smi = Chem.MolToSmiles(mol, isomericSmiles = False)
                mol = Chem.MolFromSmiles(smi)
            else:
                print("Invalid mol", seq, smi)
                continue

            proc_seqs.append(seq)
            smiles.append(smi)
            mol_list.append(mol)



        fps = self.fingerprintfn(mol_list)
        return proc_seqs, fps, smiles


    def fitness_function(self, gen, cached_dist_to_skip_calculation=None):
        """Calculates the probability of survival of each seq in generation "gen"
    
        Arguments:
            gen {list} -- sequences
            gen_n {int} -- generation number

        Returns:
            distance_av,distance_min {int} -- average and minumum distances of gen
            dist_dict, survival_dict {dict} -- {seq:distance}, {seq:probability_of_survival}
        """

        dist_dict = {}
        gen_to_calc = []

        for seq in gen:
            if cached_dist_to_skip_calculation and seq in cached_dist_to_skip_calculation:
                    dist_dict[seq] = cached_dist_to_skip_calculation[seq]
            else:
                gen_to_calc.append(seq)


        seqs, fps, smiles_l = self.calc_fingerprints(gen_to_calc)

        for i, seq in enumerate(seqs):
            map4 = fps[i]
            smiles = smiles_l[i]
            distance = self.distancefn(self.query_fp, map4)
            if distance <= self.sim_treshold:
                utils.write_results(self.results_path, smiles, seq, map4, distance)
            dist_dict[seq] = distance

        survival_dict = {}

        for k, v in dist_dict.items():
            survival_dict[k] = 1 - v 

        survival_sum = sum(survival_dict.values())
        survival_dict = {k: (v / survival_sum) for k, v in survival_dict.items()}

        distance_av = sum(dist_dict.values()) / len(dist_dict.values())
        distance_min = min(dist_dict.values())

        return distance_av, distance_min, dist_dict, survival_dict

    def who_lives(self, surv_dict):
        """Returns the sequences that will remain unchanged
        
        Returns:
            list -- chosen sequences that will live
        """

        sorted_gen = sorted(surv_dict.items(), key=lambda x: x[1], reverse=True)
        fraction = int((1 - self.gen_gap) * self.pop_size)
        if len(list(surv_dict.keys())) <= fraction:
            return list(surv_dict.keys())
        else:
            wholives = []
            if self.selec_strategy == 'Elitist':
                for element in range(fraction):
                    wholives.append(sorted_gen[element][0])
                return wholives
            elif self.selec_strategy == 'Pure':
                while len(wholives) < fraction:
                    new = np.random.choice(list(surv_dict.keys()), 1, p=list(surv_dict.values()))[0]
                    if new not in wholives:
                        wholives.append(new)
                return wholives
            else:
                if self.verbose:
                    print('not valid selection strategy, type has to be "Elitist", or "Pure"')

    def pick_parents(self, surv_dict):
        """Picks two sequences according to their survival probabilities

        Arguments:
            surv_dict {dict} -- {sequence:survival_probability}

        Returns:
            list -- parents
        """

        parents = np.random.choice(list(surv_dict.keys()), 2, p=list(surv_dict.values()))
        return parents

    def make_new_gen(self, n, surv_dict, mating_fraction, mutate):
        """Generates a new generation of n sequences with mating + 2 random sequences

        Arguments:
            n {int} -- number of structure to generate

        Returns:
            list -- new generation
        """

        new_gen = []

        for i in range(int(self.pop_size / self.rndm_newgen_fract)):
            new_gen.append(self.sequence_rng.generate())

        while len(new_gen) < n * mating_fraction:
            parents = self.pick_parents(surv_dict)
            child = utils.mating(parents)
            if self.porpouse == 'cyclic':
                child = sequence.remove_SS_cyclization(child)
            new_gen.append(child)

        # if mating fraction is less than 1.0 we fill up the gen with survivors
        while len(new_gen) < n:
            new_gen.append(np.random.choice(list(surv_dict.keys()), 1, p=list(surv_dict.values()))[0])

        if mutate:
            new_gen = mutations.mutate(new_gen, cyclic=self.porpouse == 'cyclic', mut_rate=self.mut_n, methyl=self._methyl)
        
        return new_gen


    def write_param(self):
        with open('{}/param.txt'.format(self.folder), '+w') as outFile:
            outFile.write(str(self.__dict__) + '\n')
            outFile.write('Class variables: ' + '\n')
            outFile.write('used AA: ' + str(sequence.AA) + '\n')
            outFile.write('number of point mutation: ' + str(self.mut_n) + '\n')
            outFile.write('insert branching unit rate (mutation): ' + str(self.b_insert_rate) + '\n')
            outFile.write('survival strategy: ' + str(self.selec_strategy) + '\n')
            outFile.write('fraction of new generation that is random: ' + str(self.rndm_newgen_fract) + '\n')

    def set_time_limit(self, timelimit):
        """Sets the specified timelimit. 
        the GA will stop if the timelimit is reached even if the primary condition is not reached.
        
        Arguments:
            timelimit {string} -- hours:minutes:seconds
        """

        timelimit = timelimit.split(':')
        hours = int(timelimit[0])
        minutes = int(timelimit[1])
        seconds = int(timelimit[2])
        self.timelimit_seconds = int(seconds + minutes * 60 + hours * 3600)
        if self.verbose:
            print('The GA will stop after', timelimit[0], 'hours,', timelimit[1], 'minutes, and', timelimit[2],
                  'seconds')


    def run(self):
        """Performs the genetic algorithm
 
        """
        startTime = time.time()
        
        found_identity = 0 
        
        # generation 0:
        gen = self.rndm_gen()

        if self.porpouse == 'cyclic':
            gen_cy = mutations.mutate(gen, cyclic=self.porpouse == 'cyclic',mut_rate=self.mut_n, methyl=self._methyl)
            gen = gen_cy
        if self.verbose:
            print('Generation', self.gen_n)

        # fitness function and survival probability attribution:
        distance_av, distance_min, dist_dict, surv_dict = self.fitness_function(gen, cached_dist_to_skip_calculation=None)

        if self.verbose:
            print('Average distance =', distance_av, 'Minimum distance =', distance_min)

        # progress file update (generations and their 
        # average and minimum distance from query): 
        
        utils.write_progress(self.output_path, dist_dict, self.gen_n, distance_av, distance_min)

        time_passed = int(time.time() - startTime)

        if self.verbose:
            utils.print_time(time_passed)

        # if query is found updates found identity count:
        
        if distance_min == 0:
            found_identity += 1

            # updates generation number:
        self.gen_n += 1

        # default: GA runs for ten more generation after the query is found.
        while distance_min != 0 or found_identity <= 10:

            if self.timelimit_seconds is not None and time_passed > self.timelimit_seconds:
                if self.verbose:
                    print('time limit reached')
                break

            if self.verbose:
                print('Generation', self.gen_n)

            # the sequences to be kept intact are chosen:
            survivors = self.who_lives(surv_dict)

            # n. (pop size - len(survivors)) sequences 
            # are created with crossover or mutation (cyclic):
            if self.porpouse == 'cyclic':
                new_gen = self.make_new_gen(self.pop_size - len(survivors), surv_dict, mating_fraction=0.5, mutate=True)
            else:
                new_gen = self.make_new_gen(self.pop_size - len(survivors), surv_dict, mating_fraction=1.0, mutate=False)

            # the next generation is the results of merging 
            # the survivors with the new sequences:
            gen_merg = survivors + mutations.mutate(new_gen, cyclic=self.porpouse == 'cyclic', mut_rate=self.mut_n, methyl=self._methyl, b_insert_rate=self.b_insert_rate)

            # eventual duplicates are removed:
            gen = utils.remove_duplicates(gen_merg)

            if self.verbose:
                for s in utils.random_subset(gen, 5, seed=self.seed):
                    print(f"Generated Sequence:  {sequence.reinterprete(s)}")

            # fitness function and survival 
            # probability attribution:
            distance_av, distance_min, dist_dict, surv_dict = self.fitness_function(gen, cached_dist_to_skip_calculation=dist_dict)
            if self.verbose:
                print('Average distance =', distance_av, 'Minumum distance =', distance_min)

            # progress file update (generations and their 
            # average and minimum distance from query): 
            utils.write_progress(self.output_path, dist_dict, self.gen_n, distance_av, distance_min)

            time_passed = int(time.time() - startTime)

            if self.verbose:
                utils.print_time(time_passed)

            # updates generation number (class variable):
            self.gen_n += 1

            # if query is found updates found identity count (class variable) 
            if distance_min == 0:
                found_identity += 1
