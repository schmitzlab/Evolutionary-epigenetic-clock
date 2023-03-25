'''
Nan Yao 

Oct, 2022

UGA, Athens, GA
'''



from ntpath import join
from simuOpt import setOptions
#save RAM space
#don't put it after 'import simuPOP'

setOptions(optimized=True, alleleType='binary', quiet=True)



#------------------------------------------------------------------------------------------------------#

import simuPOP as sim 
import numpy as np
import sys 

import copy
import random

from mergeGenotype import *

#------------------------------------------------------------------------------------------------------#
def simulate_selfing(popSize,lociNum, init_genotype, r, P_matrix,interval):
    pop = sim.Population(popSize,loci=lociNum,ancGen=0)
    pop.evolve(
        initOps=[sim.InitGenotype(genotype=init_genotype)],
        matingScheme=sim.SelfMating(ops=sim.Recombinator(rates=r)),
        preOps=sim.MatrixMutator(rate = P_matrix),
        gen = interval)

    output = list(pop.individual(0).genotype())
    genotype = copy.deepcopy(output)
    return genotype



def simulate_clonal(popSize,lociNum, init_genotype, r, P_matrix,interval):
    pop = sim.Population(popSize,loci=lociNum,ancGen=0)
    pop.evolve(
        initOps=[sim.InitGenotype(genotype=init_genotype)],
        matingScheme=sim.RandomSelection(ops=[sim.CloneGenoTransmitter()]),
        preOps=sim.MatrixMutator(rate = P_matrix),
        gen=interval)

    output = list(pop.individual(0).genotype())
    genotype = copy.deepcopy(output)
    return genotype




#------------------------------------------------------------------------------------------------------#



popSize = 1
#lociNum = int(1e5)


depth = 10

interval = 10


a = float(sys.argv[2])
b = float(sys.argv[3])

#recombination rate
r = 0


output_file_name_prefix = sys.argv[1]

unphased_substitution_wsp = open(output_file_name_prefix+"_substitution.csv",'a')
phased_hamming_distance_wsp = open(output_file_name_prefix+"_Hamming_distance.csv",'a')



#------------------------------------------------------------------------------------------------------#
#load initial genotype from file
#from a string to a list
#------------------------------------------------------------------------------------------------------#
init_genotype_Wsp = open(sys.argv[4])
init_genotype_str = init_genotype_Wsp.read()
init_genotype_str = init_genotype_str.strip("\n")

init_genotype_str_len = len(init_genotype_str)
init_genotype_lis = [int(init_genotype_str[i]) for i in range(init_genotype_str_len)]

print("\n\n\nload initial genotype")
print(sys.argv[4])
print("Check point")
print("Total number of loaded sites:")
print(len(init_genotype_lis))
#------------------------------------------------------------------------------------------------------#

seed=int(sys.argv[5])
sim.setRNG(seed=seed)
random.seed(seed)
print("random seed",seed)

#------------------------------------------------------------------------------------------------------#

init_genotype = [int(i) for i in init_genotype_lis]

#init_genotype_deep_copy = copy.deepcopy(init_genotype)

lociNum = int(len(init_genotype)/2)
#------------------------------------------------------------------------------------------------------#

print("population size:",popSize)
print("number of loci:",lociNum)
print("depth of simulation",depth)

print("alpha",a)
print("beta",b)

#------------------------------------------------------------------------------------------------------#
ref_genotype_Wsp = open(sys.argv[7])
ref_genotype_str = ref_genotype_Wsp.read()
ref_genotype_str = ref_genotype_str.strip("\n")

ref_genotype_str_len = len(ref_genotype_str)
ref_genotype_lis = [int(ref_genotype_str[i]) for i in range(ref_genotype_str_len)]

print("\n\n\nload reference genotype")
print(sys.argv[7])
print("Check point")
print("Total number of loaded sites:")
print(len(init_genotype_lis))

#------------------------------------------------------------------------------------------------------#

#pop = sim.Population(popSize,loci=lociNum)

P_matrix = [[1-a,a],[b,1-b]]



unphased_ref_genotype = dip2hap(ref_genotype_lis)
for i in range(depth):
    
    current_phased_temp=simulate_clonal(popSize,lociNum, init_genotype, r, P_matrix,interval)
    current_phased = copy.deepcopy(list(current_phased_temp))
    current_unphased = dip2hap(current_phased)

    #unphased
    difference_array = difference_9(unphased_ref_genotype,current_unphased)

    subsitution_out_line = ",".join([str(item) for item in difference_array]) + "\n"
    unphased_substitution_wsp.write(subsitution_out_line)


    #phased
    hamming_d,seq_len = Hamming_distance(ref_genotype_lis,current_phased)

    phased_hamming_distance_wsp.write(str(hamming_d)+","+str(seq_len)+"\n")
    
    
    #----------------------------------
    #update initial geontype
    init_genotype=current_phased







phased_genotype_out_wsp = open(sys.argv[6],'a')
phased_genotype_str = [str(item) for item in current_phased]
phased_genotype_out_line = "".join(phased_genotype_str) + "\n"
phased_genotype_out_wsp.write(phased_genotype_out_line)

print("\nphased genotype saved at\n",sys.argv[6])