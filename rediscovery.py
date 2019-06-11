from rdkit import Chem
import numpy as np
import time
import string_scoring_functions as sc
import string_crossover as co
import string_GA as ga 
import sys
from multiprocessing import Pool

Celecoxib = 'O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N'
target = Chem.MolFromSmiles(Celecoxib)
n_tries = 2
population_size = 20 
mating_pool_size = 20
generations = 20
mutation_rate = 0.5
co.average_size = target.GetNumAtoms() 
co.size_stdev = 5
co.string_type = 'smiles'
scoring_function = sc.rediscovery
max_score = 1.0
scoring_args = [target]
n_cpus = 1


file_name = sys.argv[1]

print('target', Celecoxib)
print('population_size', population_size)
print('mating_pool_size', mating_pool_size)
print('generations', generations)
print('mutation_rate', mutation_rate)
print('max_score', max_score)
print('average_size/size_stdev', co.average_size, co.size_stdev)
print('string type', co.string_type)
print('initial pool', file_name)
print('number of tries', n_tries)
print('number of CPUs', n_cpus)
print('')

results = []
size = []
t0 = time.time()
all_scores = []
generations_list = []
args = n_tries*[[population_size, file_name,scoring_function,generations,mating_pool_size,mutation_rate,scoring_args,max_score]]
with Pool(n_cpus) as pool:
    output = pool.map(ga.GA, args)

for i in range(n_tries):     
    #(scores, population) = ga.GA([population_size, file_name,scoring_function,generations,mating_pool_size,mutation_rate,scoring_args])
    (scores, population, generation) = output[i]
    all_scores.append(scores)
    print(f'{i} {scores[0]:.2f} {co.string2smiles(population[0])} {generation}')
    results.append(scores[0])
    generations_list.append(generation)
    #size.append(Chem.MolFromSmiles(sc.max_score[1]).GetNumAtoms())

t1 = time.time()
print('')
print(f'max score {max(results):.2f}, mean {np.array(results).mean():.2f} +/- {np.array(results).std():.2f}')
print(f'mean generations {np.array(generations_list).mean():.2f} +/- {np.array(generations_list).std():.2f}')
print(f'time {(t1-t0)/60.0:.2f} minutes')
#print(max(size),np.array(size).mean(),np.array(size).std())

#print(all_scores)
