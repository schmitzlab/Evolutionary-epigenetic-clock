'''
Nan Yao 

Jan, 2023

UGA, Athens, GA
'''


from ensurepip import bootstrap
from Bio import SeqIO
import sys

from math import log
from skbio import DistanceMatrix
from skbio.tree import nj

from random import choices

import numpy as np

from datetime import datetime






def S_calculator(seq1,seq2):
    #input can be str or a list of str
    length_1 =len(seq1)
    length_2 =len(seq2)

    count_0 =seq1.count("0") + seq2.count("0")
    count_1 =seq1.count("1") + seq2.count("1")

    k = count_0/count_1

    S_count = 0
    no_missing_site_count = 0

    if length_1 != length_2:
        print("lengthError")
    else:

        for i in range(length_1):
            epigenotype_1 = seq1[i]
            epigenotype_2 = seq2[i]

            if (epigenotype_1 != "?") and (epigenotype_2!="?"):
                #no missing 
                no_missing_site_count = no_missing_site_count + 1
                if epigenotype_1 != epigenotype_2:
                    #substitution
                    S_count=S_count+1
                

    S = S_count/no_missing_site_count

    return S, k

def GTR2_distance(S,k):
    term1 = (-2*k)/((1+k)**2)
    term2 = log(1+((1/term1)*S))

    return term1*term2

def GTR2_Gamma_distance(S,k,alpha):
    c = (2*k)/((1+k)**2)
    term1 = alpha*c 
    term2 = (1-(S/c))**(-1/alpha)

    return term1*(term2-1)

def GTR2_inv_Gamma_distance(S,k,alpha,p0):
    #When proportionof invariable sites, p0=0 
    #I+Gamma will be reduced to Gamma model 

    c = (2*k)/((1+k)**2)
    term1 = alpha*(1-p0)*c 
    term2 = (1-((S/c)/(1-p0)))**(-1/alpha)

    return term1*(term2-1)



#################################
def bootstrap_row_index_list_maker(seq_length):

    index_list = list(range(seq_length))
    new_index_list = choices(index_list,k=seq_length)

    return new_index_list            



#################################
#################################
#################################



current_dateTime = datetime.now()

print(current_dateTime)

seed_for_np = int(np.random.rand() * (2**32 - 1))
np.random.seed(seed=seed_for_np)

print("Random seed",seed_for_np)



fasta_file_name = sys.argv[1]
fasta_wsp = open(fasta_file_name)

GTR2_tree_file_name= sys.argv[2] + "_GTR2_tree.nwk"
GTR2_tree_wsp = open(GTR2_tree_file_name,'a')

GTR2_gamma_tree_file_name= sys.argv[2] + "_GTR2_inv_gamma_tree.nwk"
GTR2_gamma_tree_wsp = open(GTR2_gamma_tree_file_name,'a')


dis_table_file_name = sys.argv[2] + "_distance.csv"
dis_table_wsp = open(dis_table_file_name,'a')



alpha_from_input = float(sys.argv[3])
p0_from_input = float(sys.argv[5])

bootstrap_flag = sys.argv[4]



print("data file name",sys.argv[1])
print("output file prefix",sys.argv[2])
print("Alpha from input",alpha_from_input)
print("p0 from input",p0_from_input)

print("Boostrap flag",bootstrap_flag)

key_list = []
for record in SeqIO.parse(fasta_file_name, "fasta"):
    seq_name = record.id
    key_list.append(str(seq_name)) 
    #print(seq_name)

fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_wsp, "fasta"))

#the order is the same with what in the fasta file



n_taxa = len(key_list)
print("number of taxa",n_taxa)


seq_len = len(str(fasta_dict[key_list[0]].seq))
#index_projection_list = bootstrap_row_index_list_maker(seq_len)
print("sequence length",seq_len)

index_projection_array = np.random.randint(0,seq_len,seq_len)
index_projection_list = index_projection_array.tolist()

GTR2_dis_mat =[]
GTR2_gamma_dis_mat =[]
for a_key in key_list:
    print(a_key)
    
    the_seq_raw = str(fasta_dict[a_key].seq)

    if bootstrap_flag == "T":
        the_seq = [the_seq_raw[i] for i in index_projection_list]
    elif bootstrap_flag == "F":
        the_seq = the_seq_raw
    else:
        print("Error:Input bootstrap flag error!")
        print("The value of boostrap-flag must be T or F")

    #if [the_seq_raw[x] for x in range(len(the_seq_raw))] != the_seq:
    #    print("BS is done")
    

    GTR2_dis_list=[]
    GTR2_Gamma_dis_list=[]
    for another_key in key_list:

        another_seq_raw = str(fasta_dict[another_key].seq)

        if bootstrap_flag == "T":
            another_seq = [another_seq_raw[j] for j in index_projection_list]
        elif bootstrap_flag == "F":
            another_seq = another_seq_raw
        else:
            print("Error:Input bootstrap flag error!")
            print("The value of boostrap-flag must be T or F")

        this_pair_S, this_pair_k = S_calculator(the_seq,another_seq)
        #d_GTR = GTR2_distance(the_seq,another_seq)
        #d_GTR_gamma = GTR2_Gamma_distance(the_seq,another_seq,alpha=alpha_from_input)

        d_GTR = GTR2_distance(S=this_pair_S, k=this_pair_k)
        #d_GTR_gamma = GTR2_Gamma_distance(S=this_pair_S, k=this_pair_k, alpha=alpha_from_input)
        d_GTR_gamma_Inv = GTR2_inv_Gamma_distance(S=this_pair_S, k=this_pair_k, alpha=alpha_from_input, p0=p0_from_input)

        #print(a_key,another_key,d_GTR)
        #out_line= a_key + "," + another_key + "," + str(d_GTR) + "," + str(d_GTR_gamma) + "\n"
        out_line= a_key + "," + another_key + "," + str(d_GTR) + "," + str(d_GTR_gamma_Inv) + "\n"
        dis_table_wsp.write(out_line)
        
        GTR2_dis_list.append(d_GTR)
        GTR2_Gamma_dis_list.append(d_GTR_gamma_Inv)

    GTR2_dis_mat.append(GTR2_dis_list)
    GTR2_gamma_dis_mat.append(GTR2_Gamma_dis_list)




GTR2_mat = DistanceMatrix(GTR2_dis_mat,key_list)
GTR2_gamma_mat = DistanceMatrix(GTR2_gamma_dis_mat,key_list)

GTR2_tree_object = nj(GTR2_mat)
GTR2_gamma_tree_object = nj(GTR2_gamma_mat)

GTR2_tree_object.write(GTR2_tree_wsp,'newick')
print("Tree based on GTR2 distance")
print(GTR2_tree_object.ascii_art())

GTR2_gamma_tree_object.write(GTR2_gamma_tree_wsp,'newick')
print("Tree based on GTR2+I+G distance")
print(GTR2_gamma_tree_object.ascii_art())

print("Compare tree topologies")
print("Robinson and Foulds symmetric difference:")
print(GTR2_gamma_tree_object.compare_rfd(GTR2_tree_object))

print("\n\n------------------\n\n")