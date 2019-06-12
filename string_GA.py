'''
Written by Jan H. Jensen 2018. 
Many subsequent changes inspired by https://github.com/BenevolentAI/guacamol_baselines/tree/master/graph_ga
'''

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdmolops

from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

import numpy as np
import random
import time
import sys

import string_crossover as co
import string_mutate as mu
import string_scoring_functions as sc

def read_file(file_name):
  smiles_list = []
  with open(file_name,'r') as file:
    for smiles in file:
      smiles = smiles.replace('\n', '')
      smiles_list.append(smiles)

  return smiles_list

def make_initial_population(population_size,file_name):
  smiles_list = read_file(file_name)
  population = []
  for _ in range(population_size):
    smiles = random.choice(smiles_list)
    mol = Chem.MolFromSmiles(smiles)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
    string = co.smiles2string(smiles)
    if string:
      population.append(string)
    
  return population

def calculate_normalized_fitness(scores):
  sum_scores = sum(scores)
  normalized_fitness = [score/sum_scores for score in scores]

  return normalized_fitness

def make_mating_pool(population,fitness,mating_pool_size):
  mating_pool = []
  for i in range(mating_pool_size):
  	mating_pool.append(np.random.choice(population, p=fitness))

  return mating_pool
 

def reproduce(mating_pool,population_size,mutation_rate):
  new_population = []
  while len(new_population) < population_size:
    parent_A = random.choice(mating_pool)
    parent_B = random.choice(mating_pool)
    new_child = co.crossover(parent_A,parent_B)
    if new_child != None:
	    new_child = mu.mutate(new_child,mutation_rate)
	    if new_child != None:
	    	new_population.append(new_child)

  return new_population

def sanitize(population,scores,population_size,prune_population):
    if prune_population:
      smiles_list = []
      population_tuples = []
      for score, string in zip(scores,population):
          canonical_smiles = Chem.MolToSmiles(co.string2mol(string))
          if canonical_smiles not in smiles_list:
              smiles_list.append(canonical_smiles)
              population_tuples.append((score, string))
    else:
      population_tuples = list(zip(scores,population))


    population_tuples = sorted(population_tuples, key=lambda x: x[0], reverse=True)[:population_size]
    new_population = [t[1] for t in population_tuples]
    new_scores = [t[0] for t in population_tuples]

    return new_population, new_scores

def GA(args):
  population_size, file_name, scoring_function, generations, mating_pool_size, mutation_rate, \
  scoring_args, max_score, prune_population = args
 
  population = make_initial_population(population_size,file_name)
  scores = sc.calculate_scores(population,scoring_function,scoring_args)
  fitness = calculate_normalized_fitness(scores)

  for generation in range(generations):
    mating_pool = make_mating_pool(population,fitness,mating_pool_size)
    new_population = reproduce(mating_pool,population_size,mutation_rate)
    new_scores = sc.calculate_scores(new_population,scoring_function,scoring_args)
    population, scores = sanitize(population+new_population, scores+new_scores, population_size,prune_population)  
    fitness = calculate_normalized_fitness(scores)
    if scores[0] >= max_score:
      break

  return (scores, population, generation+1)


if __name__ == "__main__":
    Celecoxib = 'O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N'
    target = Chem.MolFromSmiles(Celecoxib)
    population_size = 20 
    mating_pool_size = 20
    generations = 20
    mutation_rate = 0.01
    co.average_size = 39.15
    co.size_stdev = 3.50
    co.string_type = 'selfies'
    scoring_function = sc.rediscovery #sc.logP_score
    #scoring_function = sc.logP_score
    max_score = 1.0 # 9999.
    prune_population = True
    scoring_args = [target]

    file_name = 'ZINC_first_1000.smi'

    (scores, population, generation) = GA([population_size, file_name, scoring_function, generations,
                                           mating_pool_size, mutation_rate, scoring_args, max_score,
                                           prune_population])
    
    print('done')
