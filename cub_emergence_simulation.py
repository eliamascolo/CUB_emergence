# -*- coding: utf-8 -*-
"""
Simulation of emerging Codon Usage Bias.

"""

import random
import numpy as np
import json
import copy


class Genome():
    
    def __init__(self, config_dict):
        
        # Set parameters from config file
        self.n_codon_types = config_dict['number_of_synonymous_codons']
        self.n_aa_occurences = config_dict['number_of_aa_occurences_in_genome']
        self.wait_time_cost_param = config_dict['wait_time_cost_param']
        self.max_trna_freq_shift = config_dict['max_trna_freq_shift']
        self.trna_freq_lower_bound = config_dict['trna_freq_lower_bound']
        self.trna_freq_upper_bound = (
            1 - (self.trna_freq_lower_bound * (self.n_codon_types - 1))
        )
        self.max_error_trna_freq_sum = config_dict['max_error_trna_freq_sum']
        self.mutate_trna_probability = config_dict['mutate_trna_probability']
        
        # Codon counts
        self.codon_counts = []
        
        # tRNA aboundances
        self.trna_abundances = []
        
        # Fitness
        self.fitness = None
    
    def assign_codons_evenly(self):
        '''
        Synonymous codons are represented as evenly as possible.
        n_aa_occurences:
            the number of times the a.a. is present
        n_codon_types:
            the number of synonymous codons by which the a.a. can be encoded
        The function sets the `codon_counts` attribute.
        self.codon_counts:
            A list containing the number of instances for each codon type.
            Therefore,
            len(self.codon_counts) == self.n_codon_types
            and
            sum(self.codon_counts) == self.n_aa_occurences
        
        '''
        quotient = self.n_aa_occurences // self.n_codon_types
        remainder = self.n_aa_occurences % self.n_codon_types
        
        small_counts = [quotient] * (self.n_codon_types - remainder)
        larger_counts = [quotient + 1] * remainder
        counts = small_counts + larger_counts
        
        # Set codon counts
        self.codon_counts = counts
    
    def assign_tRNAs_abundances_evenly(self):
        '''
        Sets the abundances of tRNAs recognizing synonymous codons to equal
        values. The function sets the `trna_abundances` attribute.
        self.trna_abundances:
            A list containing the abundance value for each tRNA type.
        '''
        n_tRNAs = self.n_codon_types
        
        # Each tRNA will have a frequency equal to freq
        freq = 1 / n_tRNAs
        
        # Set tRNA aboundances
        self.trna_abundances = [freq] * n_tRNAs
    
    
    def tRNA_abundance_is_within_bounds(self, trna_abundance):
        
        return (trna_abundance > self.trna_freq_lower_bound and
                trna_abundance < self.trna_freq_upper_bound)
    
    def mutate(self):
        
        if random.random() < self.mutate_trna_probability:
            self.mutate_relative_tRNA_abundance()
        else:
            self.mutate_codon()
    
    def mutate_codon(self):
        '''
        Produces a synonymous mutation, changing one random condon into a
        synonymous codon.
        '''
        donor, acceptor = random.sample(range(self.n_codon_types), 2)
        
        if self.codon_counts[donor] > 0:
            
            # One codon mutates to a synonymous codon 
            self.codon_counts[donor] -= 1
            self.codon_counts[acceptor] += 1
    
    def mutate_relative_tRNA_abundance(self):
        
        random_tRNA_type = random.randint(0, self.n_codon_types - 1)
        current_abundance = self.trna_abundances[random_tRNA_type]
        
        s = self.max_trna_freq_shift
        
        shift = random.uniform(-s, s)
        
        if current_abundance + shift < self.trna_freq_lower_bound:
            new_abundance = self.trna_freq_lower_bound
        elif current_abundance + shift > self.trna_freq_upper_bound:
            new_abundance = self.trna_freq_upper_bound
        else:
            new_abundance = current_abundance + shift
        
        # Set new relative abundance for the selected tRNA type
        self.trna_abundances[random_tRNA_type] = new_abundance
        
        # Adjust the other relative abundances
        sum_of_the_other_abundances = 1 - current_abundance
        
        scaling_factor = 1 - (shift / sum_of_the_other_abundances)
        
        for i in range(self.n_codon_types):
            if i != random_tRNA_type:
                
                if self.tRNA_abundance_is_within_bounds(
                        self.trna_abundances[i] * scaling_factor
                ):
                    self.trna_abundances[i] *= scaling_factor
                
                else:
                    print("tRNA freq went out of bounds. Debug this please")
        
        error = sum(self.trna_abundances) - 1
        
        if abs(error) > self.max_error_trna_freq_sum:
            
            for i in range(self.n_codon_types):
                # try to remove the error
                trna_abundance = self.trna_abundances[i]
                if (not trna_abundance - error < self.trna_freq_lower_bound and
                    not trna_abundance - error > self.trna_freq_upper_bound):
                    
                    self.trna_abundances[i] -= error
                    break
                
                else:
                    continue
    
    def evaluate_fitness(self):
        
        total_cost = 0
        
        for i in range(self.n_codon_types):
            codon_class_cost = self.get_cost_for_codon_class(i)
            total_cost += codon_class_cost
        
        self.fitness = - total_cost
        
    def get_cost_for_codon_class(self, codon_type : int):
        
        # Abundance of that codon class
        codon_abundance = self.codon_counts[codon_type]
        
        # Relative abundance (frequency) of the cognate tRNA
        freq_of_cognate_trna = self.trna_abundances[codon_type]
        
        expected_cost = self.get_expected_wait_time_cost(freq_of_cognate_trna)
        
        codon_class_cost = expected_cost * codon_abundance
        
        return codon_class_cost
    
    # def get_expected_wait_time_cost(self, tRNA_freq):
    #     expected_wait_time = 1 / (-np.log(1 - tRNA_freq))
    #     expected_cost =  self.wait_time_cost_param * expected_wait_time
    #     return expected_cost
    
    def get_expected_wait_time_cost(self, tRNA_freq):
        
        return 1 / tRNA_freq
    
    def replicate(self):
        
        clone = copy.deepcopy(self)
        
        return clone



def read_json_file(filename):
    with open(filename) as json_content:
        return json.load(json_content)



config_dictionary = read_json_file('config.json')


# TEST SPACE

new_genome = Genome(config_dictionary)
new_genome.assign_codons_evenly()
new_genome.assign_tRNAs_abundances_evenly()

new_genome2 = new_genome.replicate()







































