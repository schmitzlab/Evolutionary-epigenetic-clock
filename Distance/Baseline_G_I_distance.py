'''
Nan Yao 

Jan, 2023

UGA, Athens, GA
'''



import scipy.optimize as optimize
import numpy as np 

from scipy.linalg import expm
from math import log

import sys


def loadsubs(fileName):

    wsp = open(fileName)

    lines = wsp.readlines()

    subMat = []

    for i in range(len(lines)):


        line = lines[i]
        line = line.strip('\n')

        eachIndividual = line.split(',')

        eachIndividual = [int(item) for item in eachIndividual]

        #length of eachIndivdual is 9
        subMat.append(eachIndividual)

    print("substitution file has been loaded")
    print("file name:",fileName)
    print("number of lines:",len(lines))
    print("length of each line:",len(subMat[0]))

    return subMat







def tranProbFunc_baseline(para):

    a_t, k = para
    
    Qt = np.asarray([[-2*a_t, 2*a_t, 0],
                     [k*a_t,(-1-k)*a_t,a_t],
                     [0, 2*k*a_t, -2*k*a_t]])

    P_t = expm(Qt)
    out = P_t.reshape((9,))

    return out.tolist() 


def gamma_inv_MGF(lambda_k_t,p0,alpha):
    
    beta = (1-p0)*alpha
    
    gamma_mgf_value = (1-lambda_k_t/beta)**(-alpha)

    mgf_value = p0 + (1-p0)*gamma_mgf_value

    return mgf_value

#don't use this one
def tranProbFunc_Gamma_Inv_baseline(para, p0, alpha):
    #when p0 and alpha are 0 
    #the function will be reduced to baseline uniform P(t) function

    a_t, k = para
    
    c = (-1/k)
    Z = 1/((1+k)**2)

    U = [[1,c**2,c],
         [1,c,(k-1)/(2*k)],
         [1,1,1]]
    
    U_inv = [[(k**2)*Z, (2*k)*Z, Z],
             [(k**2)*Z, -2*(k**2)*Z, (k**2)*Z],
             [-2*(k**2)*Z, (2*k)*(k-1)*Z, (2*k)*Z]]

    lambda_vec = [0, -2*(1+k)*a_t, -1*(1+k)*a_t ]

    P_mat = []

    for state_i in range(3):

        P_i_list = []

        for state_j in range(3):

            P_ij = 0

            for state_k in range(3):

                u_ik = (U[state_i])[state_k]
                u_inv_kj = (U_inv[state_k])[state_j]

                mgf_value = gamma_inv_MGF(lambda_k_t=lambda_vec[state_k],p0=p0,alpha=alpha)

                _temp = u_ik*u_inv_kj*mgf_value

                P_ij = P_ij + _temp

            P_i_list.append(P_ij)

        P_mat.append(P_i_list)

    P_t = np.array(P_mat)

    out = P_t.reshape((9,))

    return out.tolist()



###########
def tranProbFunc_Gamma_Inv_baseline_v2(para, p0, alpha):
    #when p0 and alpha are 0 
    #the function will be reduced to baseline uniform P(t) function

    a_t, k = para
    
    c = (-1/k)
    Z = 1/((1+k)**2)

    U = [[1,c**2,c],
         [1,c,(k-1)/(2*k)],
         [1,1,1]]

    U = np.array(U)
    
    U_inv = [[(k**2)*Z, (2*k)*Z, Z],
             [(k**2)*Z, -2*(k**2)*Z, (k**2)*Z],
             [-2*(k**2)*Z, (2*k)*(k-1)*Z, (2*k)*Z]]

    U_inv = np.array(U_inv)

    lambda_vec = [0, -2*(1+k)*a_t, -1*(1+k)*a_t ]

    mgf_vec = [gamma_inv_MGF(lambda_k_t=L,p0=p0,alpha=alpha) for L in lambda_vec]

    D = np.zeros((3, 3))
    np.fill_diagonal(D, mgf_vec)

    term_1 = np.matmul(U,D)
    P_t = np.matmul(term_1,U_inv)

    out = P_t.reshape((9,))

    return out.tolist()



###########





def totalMiu(para):

    a_t, k = para

    #here was a bug

    coef = (4*k)/(k+1)
    return coef * a_t






def negLogLL(data,p0,alpha):


    def logLLWithData(para):


        #numericPt = transProbFunc(para)
        numericPt = tranProbFunc_Gamma_Inv_baseline_v2(para=para,p0=p0,alpha=alpha)

        a_t, k = para
        
        pi_1 = (k**2)/((1+k)**2)
        pi_2 = (2*k)/((1+k)**2)
        pi_3 = 1 - pi_1 -pi_2

        #piList = [pi_1,pi_2,pi_3, pi_1,pi_2,pi_3, pi_1,pi_2,pi_3]

        #Fixed on Sep 12 
        piList = [pi_1,pi_1,pi_1, pi_2,pi_2,pi_2, pi_3,pi_3,pi_3]


        #To keep it safe for log operation
        for i in range(9):
            if numericPt[i] <=0:
                numericPt[i] = 10**-10

            if piList[i] <=0:
                piList[i] = 10**-10


        piPtList = [log(piList[i]*numericPt[i]) for i in range(9)]

        logLikelihoodList = [data[i]*piPtList[i] for i in range(9)]

        logLLValue = sum(logLikelihoodList)

        return -logLLValue

    return logLLWithData



subsFileName = sys.argv[1]
outputFileName = sys.argv[2]

p0 = float(sys.argv[3])
alpha = float(sys.argv[4])

subs = loadsubs(subsFileName)

dataWsp = open(subsFileName)
outputWsp = open(outputFileName,'a')

print("Input substitution table",subsFileName)
print("Output distance table",outputFileName)

print("Input p0=",p0)
print("Input alpha=",alpha)


#a_t, k
initPara = [0.001,1]

#a_t, k
bounds = [(10**-7,3),(10**-7,10**8)]

miuPerYear = 1


for i in range(len(subs)):

    print("current index is",i)

    data = subs[i]

    #GTR model with current pairwise substitution
    #targetFunc = negLogLL(data, tranProbFunc_baseline)
    targetFunc = negLogLL(data, p0, alpha)

    result = None

    result = optimize.minimize(targetFunc, initPara, bounds=bounds)

    out = None 

    if result.success:
        out = result.x
    
        #at_est, k1_est, k2_est, pi1_est, pi2_est = out

        llhValue = - targetFunc(out)

        print(out)
        
        print("llh =",llhValue)

        print("brach length = ", totalMiu(out))

        print("divergence time" , totalMiu(out)/miuPerYear)

        outLine = ""
        outLine = str(llhValue) + "," + str(totalMiu(out)) + "," + str(totalMiu(out)/miuPerYear) + ","

        estStr = ",".join([str(item) for item in out])

        outLine = outLine + estStr + "\n"
        outputWsp.write(outLine)

    else:
        #raise ValueError(result.message)
        
        #this value error will be detected when building tree
        outputWsp.write("9999,9999,9999\n")
        print("ValueError!")
        print(result.message)


    print("\n\n\n")






