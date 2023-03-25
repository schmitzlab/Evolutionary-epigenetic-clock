'''
Nan Yao 

Oct, 2022

UGA, Athens, GA
'''

from simuOpt import setOptions
#save RAM space
#don't put it after 'import simuPOP'

setOptions(optimized=True, alleleType='short', quiet=True)




#------------------------------------------------------------------------------------------------------#

import simuPOP as sim 
import numpy as np
import sys 

import copy
import random

from mergeGenotype import *

#------------------------------------------------------------------------------------------------------#
def simulate(popSize,lociNum, init_genotype, r, P_matrix,interval):
    pop = sim.Population(popSize,loci=lociNum,ancGen=0)
    pop.evolve(
        initOps=[sim.InitGenotype(genotype=init_genotype)],
        matingScheme=sim.SelfMating(ops=sim.Recombinator(rates=r)),
        preOps=sim.MatrixMutator(rate = P_matrix),
        gen = interval)

    output = list(pop.individual(0).genotype())
    genotype = copy.deepcopy(output)
    return genotype


#------------------------------------------------------------------------------------------------------#



popSize = 1
#lociNum = int(1e5)


depth = 10

interval = 10

output_file_name_prefix = sys.argv[1]

a = float(sys.argv[2])
b = float(sys.argv[3])

#recombination rate
r = 0
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

#the a and be here are actually alpha and beta 
di = 1-a-b-b
#k2p
P_matrix = [[di,a,b,b],
            [a,di,b,b],
            [b,b,di,a],
            [b,b,a,di]]

homo_list=[]
S2_list=[]

unphased_ref_genotype = dip2hap_SNP(ref_genotype_lis)
for i in range(depth):
    '''
    pop = sim.Population(popSize,loci=lociNum)
    pop.evolve(
        initOps=[sim.InitGenotype(genotype=init_genotype)],
        matingScheme=sim.SelfMating(ops=sim.Recombinator(rates=r)),
        preOps=sim.MatrixMutator(rate = P_matrix),
        gen = interval)

    current_phased_temp = pop.individual(0).genotype()
    '''
    current_phased_temp=simulate(popSize,lociNum, init_genotype, r, P_matrix,interval)
    current_phased = copy.deepcopy(list(current_phased_temp))
    current_unphased = dip2hap_SNP(current_phased)
    

    homo_count,S2_count = homo_P_distance_SNP(unphased_ref_genotype,current_unphased)
    print("homo",homo_count,"S2",S2_count)

    init_genotype=current_phased

    homo_list.append(homo_count)
    S2_list.append(S2_count)
    #del pop

S2_wsp = open(output_file_name_prefix+"_S2_counts.csv",'a')
S2_list_str = [str(item) for item in S2_list]
S2_out_line = ",".join(S2_list_str) + "\n"
S2_wsp.write(S2_out_line)

homo_count_wsp = open(output_file_name_prefix+"_homo_counts.csv",'a')
homo_list_str = [str(item) for item in homo_list]
homo_out_line = ",".join(homo_list_str) + "\n"
homo_count_wsp.write(homo_out_line)

phased_genotype_out_wsp = open(sys.argv[6],'a')
phased_genotype_str = [str(item) for item in current_phased]
phased_genotype_out_line = "".join(phased_genotype_str) + "\n"
phased_genotype_out_wsp.write(phased_genotype_out_line)

print("\nphased genotype saved at\n",sys.argv[6])

if phased_genotype_out_line == init_genotype_str:
    print("Error, updating failed")
else:
    print("input initial genotype str length",len(init_genotype_str))
    print("output current genotype str length",len(phased_genotype_out_line))
    print("\nDone-------------------------------------\n")