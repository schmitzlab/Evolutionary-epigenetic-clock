'''
Nan Yao 

Oct, 2022

UGA, Athens, GA
'''


from simuOpt import setOptions
#save RAM space
#don't put it after 'import simuPOP'

setOptions(optimized=True, alleleType='binary', quiet=True)


#------------------------------------------------------------------------------------------------------#

import simuPOP as sim 
import numpy as np
import sys 

from mergeGenotype import *

popSize = 1
lociNum = int(1e5)

depth = 500

a = float(sys.argv[2])
b = float(sys.argv[3])

#recombination rate
r = 0
#------------------------------------------------------------------------------------------------------#
#load initial genotype from file
#------------------------------------------------------------------------------------------------------#
init_genotype_Wsp = open(sys.argv[4])
init_genotype_str = init_genotype_Wsp.read()
init_genotype_str = init_genotype_str.strip("\n")

init_genotype_str_len = len(init_genotype_str)
init_genotype_lis = [int(init_genotype_str[i]) for i in range(init_genotype_str_len)]

print("Check point")
print("Total number of loaded sites:")
print(len(init_genotype_lis))

init_genotype = [int(i) for i in init_genotype_lis]
#------------------------------------------------------------------------------------------------------#

print("population size:",popSize)
print("number of loci:",lociNum)
print("depth of simulation",depth)

print("alpha",a)
print("beta",b)

#------------------------------------------------------------------------------------------------------#


outputWsp = open(sys.argv[1],'a')

pop = sim.Population(popSize,loci=lociNum,ancGen=-1)

#the a and be here are actually alpha and beta 
P_matrix = [[1-a,a],[b,1-b]]





#initOps=[sim.InitGenotype(freq=[(b/(a+b)),(a/(a+b))])],
pop.evolve(
    initOps=[sim.InitGenotype(genotype=init_genotype)],
    matingScheme=sim.SelfMating(ops=sim.Recombinator(rates=r)),
    preOps=sim.MatrixMutator(rate = P_matrix),
    gen = depth)


print('simulation done!')

#sim.dump(pop)

print("saving genotypes")


#generation by generation 
#the first row is current sequence
for i in range(depth+1):

    idx = i
    #print(idx,"generations ago")
    
    pop.useAncestralGen(idx)
    siteStates = pop.individual(0).genotype()

    #unphased
    genotypes = dip2hap(siteStates)
    
    #phased
    #genotypes = list(siteStates)
    genotypesStr = ",".join(genotypes) + "\n"

    outputWsp.write(genotypesStr)
    

    #print(siteStates[0:lociNum])
    #print(siteStates[lociNum:lociNum*2])
    #print(genotypes)

    #print(genotypesStr)


    print("")


print("simulation results have been saved to",sys.argv[1])



