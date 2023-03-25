'''
Nan Yao 

Oct, 2022

UGA, Athens, GA
'''

def loadSimuMethylomes(fileName):

    wsp = open(fileName)

    lines = wsp.readlines()

    popMat = []

    for i in range(len(lines)):


        line = lines[i]
        line = line.strip('\n')

        eachIndividual = line.split(',')

        eachIndividual = [int(item) for item in eachIndividual]

        popMat.append(eachIndividual)

    print("Methylome file loading finished \nfile name:", fileName)
    print("popualtion size is", len(popMat))
    print("length of sequence is", len(popMat[0]))

    return popMat



def dip2hap(rawSeq):

    seqLen = len(rawSeq)

    rawSeq = [str(site) for site in rawSeq]

    #print("length of raw sequence",seqLen)
    realLen = int(seqLen/2)

    epigenotypes = [] 
    for i in range(realLen):

        allele_idx_1 = i 
        allele_idx_2 = i + realLen

        allele_1 = rawSeq[allele_idx_1]
        allele_2 = rawSeq[allele_idx_2] 

        this_epigenotype = set([allele_1,allele_2])

        if this_epigenotype == set(['0','0']):

            temp = '0'

        elif this_epigenotype == set(['1','0']):

            temp = '1'

        elif this_epigenotype == set(['1','1']):

            temp = '2'

        else:

            print("allele type error!")

            temp = '-1'

        epigenotypes.append(temp)


    return epigenotypes

