#!/usr/bin/env python
import os
import sys
import numpy as np
import dpFuncs_sw2 as dpF
from Bio import SeqIO

def tb_prefix(query,prefix,P,P_p,iPre):
    """
    The objective of this function is to traceback the Prefix score matrix and return alignement sequences for prefix in reference and query
 
    Paramters
    ---------
    query = array of int, shape(n_bases_q,)
            The query sequence information

    prefix = array of int, shape(n_bases_p,)
             The prefix sequence information for the reference

    P = array of int, shape(n_bases_q,n_bases_p)
        The prefix score matrix

    P-p = array of int, shape(n_bases_q,n_bases_p)
          The prefix pointer matrix

    iPre = int
           row index where traceback jumps from TR repeat matrix R into prefix score matrix P

    Returns
    -------
    q_pre_Aligned = list
                    A list ['prefix sequence in ref','prefix sequence in query'] 

    """
    pre = False
    j = np.size(prefix)
    i = iPre
    qSeq = map(chr,query)
    preSeq = map(chr,prefix)
    qList = []
    preList= []
    jLessThanZero = False
    while(j > 0):
      if j > 0: 
         if (P[i][j] > 0):
            if (P_p[i][j] == 0):
               qList.append(qSeq[i-1])
               preList.append(preSeq[j-1])
               i += -1
               j += -1
            if (P_p[i][j] == 1):
               qList.append(qSeq[i-1])
               preList.append("-")
               i += -1
            if (P_p[i][j] == 2):
               qList.append("-")
               preList.append(preSeq[j-1])
               j += -1
         else:
            jLessThanZero = True
            break
      else:
         jLessThanZero = True
         break
    if jLessThanZero == False:
       preList.reverse()
       qList.reverse()
       q_pre_Aligned = []
       q_pre_Aligned.append(''.join(preList))
       q_pre_Aligned.append(''.join(qList))
    else:
       q_pre_Aligned = ['','']
    return q_pre_Aligned

def tb_suffix(query,suffix,S,S_p):
    """
    The objective of this function is to traceback the Suffix score matrix.
 
    Paramters
    ---------
    query = array of int, shape(n_bases_q,)
            The query sequence information

    suffix = array of int, shape(n_bases_s,)
             The suffix sequence information for the reference

    S = array of int, shape(n_bases_q,n_bases_s)
        The suffix score matrix

    S_p = array of int, shape(n_bases_q,n_bases_s)
          The suffix pointer matrix

    Returns
    -------
    index_i_tR = int
                 The index in query where suffix traceback enters the repeat matrix R

    bestscoresuffix = int
                      The maximum smith-waterman score

    q_suf_Aligned = list
                    A list ['suffix sequence in ref','suffix sequence in query'] 

    """

    temp = 0
    i = 0
    j = np.size(suffix)
    nq = np.size(query)
    qSeq = map(chr,query)
    sufSeq = map(chr,suffix)
    qList = []
    sufList= []
    for k in range(nq+1): # Finds the maximal value for the suffix to begin tb
        if temp < S[k][j]:
           temp = S[k][j]
           i = k
    bestscoresuffix = temp
    pre = False
    index_i_tR = 0
    jLessThanZero = False
    while (not(pre)):
      if j > 0:
         if (S[i][j] > 0):
            if (S_p[i][j] == 1):
               qList.append(qSeq[i-1])
               sufList.append("-")
               i += -1
            elif (S_p[i][j] == 2):
               qList.append("-")
               sufList.append(sufSeq[j-1])
               j += -1
            elif (S_p[i][j] == 0):
               qList.append(qSeq[i-1])
               sufList.append(sufSeq[j-1])
               i += -1
               j += -1
            elif (S_p[i][j] == 3):
               qList.append(qSeq[i-1])
               sufList.append(sufSeq[j-1])
               j += -1
               i += -1
               pre = True
            index_i_tR = i
         else :
            index_i_tR = 0
            jLessThanZero = True
            break
      else:
         index_i_tR = 0
         jLessThanZero = True
         break
    if jLessThanZero == False:
       sufList.reverse()
       qList.reverse()
       q_suf_Aligned = []
       q_suf_Aligned.append(''.join(sufList))
       q_suf_Aligned.append(''.join(qList))
       #print "q_suf_Aligned:\n",q_suf_Aligned[0],"\n",q_suf_Aligned[1] 
    else:
       q_suf_Aligned = ['','']
    return index_i_tR, bestscoresuffix,q_suf_Aligned

def finalTrace(qArray,tR,T,T_p,tR_rowMax,index_i_tR):
    """
    The objective of this function is to traceback the repeat score matrix and identify the TR sequence in query
 
    Paramters
    ---------
    qArray = array of int, shape(n_bases_q,)
            The query sequence information

    tR = string 
         The TR element  sequence information for the reference

    T = array of int, shape(n_bases_q + 1,n_bases_tR + 1)
        The repeat score matrix

    T_p = array of int, shape(n_bases_q + 1,n_bases_tR + 1)
          The repeat pointer matrix

    tR_rowMax = numpy array of int, shape(n_bases_q + 1, 2)
                stores the information on the maximum smith-waterman score for each row

    index_i_tR = int
                 The index from where the traceback starts in the repeat score matrix, T

    Returns
    -------
    number_tR = float
                Estimated number of TR elements in the query based on dynamic programming output

    repeatAligned = list
                    A list ['tandem repeat sequence in ref','tandem repeat sequence in query'] 

    scoretRentry = int
                   Score when entered into the repeat matrix, R

    scoretRexit = int
                  Score when exited from the repeat matrix into prefix score matrix

    naiveTR = float
              naive estimation of number of repeat elements based on estimated prefix and suffix boundaries

    """

    i = index_i_tR #This is from where the traceback starts
    #print "i index in TR: ", i
    if (i > 0):
      j = tR_rowMax[i][1]
      scoretRentry = T[i][j]
      t = len(tR)
      pre = False
      trScoreNeg = False
      query = map(chr,qArray)
      qList = []
      tList = []
      while (not (pre)):
        if j > 0:
          if (T[i][j] > 0):
            if (T_p[i][j] == 0):
               qList.append(query[i-1])
               tList.append("-")
               i += - 1
            elif (T_p[i][j] == 1):
               qList.append("-")
               tList.append(tR[j-1])
               j += - 1
            elif (T_p[i][j] == 2):
               qList.append(query[i-1])
               tList.append(tR[j-1])
               i += - 1
               j += - 1
            elif (T_p[i][j] == 3):
               qList.append(query[i-1])
               tList.append(tR[j-1])
               for k in range(0,j-1): #open up gaps on query
                   qList.append("-")
                   tList.append(tR[k])
               i += - 1
               j = t
            elif (T_p[i][j] == 4):
               qList.append(query[i-1])
               tList.append(tR[j-1])
               i += - 1
               pre = True
               trSinQ = i
               scoretRexit = T[i+1][j]
          else:
              trScoreNeg = True
              print "Tandem repeat score <= 0"
              break
        else:
            trScoreNeg = True
            break
      #print "trScoreNeg: ",trScoreNeg
      if trScoreNeg == False: #Repeat matrix score was never less than zero, implies all is well
         #print "tandem repeat score > 0"
         tList.reverse()
         qList.reverse()
         repeatAligned = []
         repeatAligned.append(''.join(tList))
         repeatAligned.append(''.join(qList))

         tempTList = []
         for i in range(0, len(tList)):
            if tList[i] == 'A' or tList[i] == 'C' or tList[i] == 'G' or tList[i] == 'T':
               tempTList.append(tList[i])

         tRLength = len(tempTList)
         str_tR = "".join(tR)
         str_tempTList = "".join(tempTList)
         if (t < tRLength):
              if t >= 1:
                 number_tR = (tRLength)/float(t) #This is the estimated number of tanderm repeat elements from the dynamic programming method
                 naiveTR = index_i_tR - trSinQ #naive number of tandem repeat elements based on boundary estimation
              else:
                 number_tR = -30000.0
                 naiveTR = -30000.0
                 repeatAligned = ['','']
                 scoretRentry = 0
                 scoretRexit = 0
         else:
           number_tR = -20000.0
           naiveTR = -20000.0
           repeatAligned = ['','']
           scoretRentry = 0
           scoretRexit = 0
      else:
          number_tR = -40000.0
          naiveTR = -40000.0
          repeatAligned = ['','']
          scoretRentry = 0
          scoretRexit = 0
    else:
      number_tR = -10000.0
      naiveTR = -10000.0
      repeatAligned = ['','']
      scoretRentry = 0
      scoretRexit = 0
    return number_tR, repeatAligned, scoretRentry, scoretRexit, naiveTR
 
def alignRegions(query,pre_suf_tR):
    """
    The objective of this function is to identify the TR boundaries within the reads (Section 2.1)

    Parameters
    ---------- 
    query = string 
           The read sequence

    pre_suf_tR = List
                 A list of strings: [<prefix_seq>, <suf_seq>, <tr_seq>]

    Returns
    -------  
    #return number_tR, repeatAligned, Prefix_BestScore, Suffix_BestScore, tR_maxScore, tR_exitScore, alignSufxQ, alignPrfxQ, naiveTR
    number_tR = float
                number of tandem repeats in query as estimated by the 3-stage dynamic programming algorithm

    repeatAligned = list
                    A list of strings: [<TR align sequence in ref>, <TR align sequence in query>]

    Prefix_BestScore = int
                       Smith-waterman score of prefix score matrix, P (Prefix)

    Suffix_BestScore = int
                       Smith-waterman score in suffix score matrix, S (Prefix + TR + Suffix)

    tR_maxScore = int
                  Smith-waterman score in TR matrix, R (Prefix + TR)

    tR_exitScore = int
                   Smith-waterman score in TR matrix where alignment jumps to prefix matrix, P

    alignSufxQ = list
                 A list of Aligned suffix sequence in read ['suffix align sequence in reference','suffix align sequence in query']

    alignPrfxQ = list
                 A list of aligned prefix sequence in read ['prefix align sequence in reference','prefix align sequence in query']

    naiveTR = float
              number of tandem repeats in query as estimated by naive estimation of prefix and suffix boundaries in read
 
    Notes
    -----
    This function implements section 2.1 from the paper. For speed ups the main dynamic programming parts have been implemented in cython.
    """
    temp_trfSeq = []
    prefix = []
    prefix = list(pre_suf_tR[0])
    len_prefix = len(prefix)
    pArray = np.arange(len_prefix,dtype=np.int_)
    tempP = map(ord,prefix)
    tempP2 = np.asarray(tempP)
    pArray = tempP2 #A numpy array of the PREFIX
#    print "pArray: ", pArray

    suffix = []
    suffix = list(pre_suf_tR[1])
    len_suffix = len(suffix)
    sArray = np.arange(len_suffix,dtype=np.int_)
    tempS = map(ord,suffix)
    tempS2 = np.asarray(tempS)
    sArray = tempS2 #A numpy array of the SUFFIX
#    print "sArray: ", sArray

    tR = []
    tR = list(pre_suf_tR[2])
    len_tR = len(tR)
    tRArray = np.arange(len_tR,dtype=np.int_)
    temptR = map(ord,tR)
    temptR2 = np.asarray(temptR)
    tRArray = temptR2 #A numpy array of the TR element
#    print "tRArray: ", tRArray

    #initialization for prefix
    preScore = np.zeros((len(query)+1,len_prefix+1),dtype=np.int_)
    # Added in initialization for prefix
    preMaxScore = np.zeros((len(query)+1,len_prefix+1),dtype=np.int_)
    preMaxScore.fill(-1)
    for j in xrange (1, len_prefix+1):
        preScore[0][j] = -5*j

    Prefix_BestScore = 0
    Prefix_BestScoreIndex = 0
    Prefix_BestScore, Prefix_BestScoreIndex = dpF.sw_prefix(query,pArray,preScore,preMaxScore) #Calling sw_prefix function (implemented in Cython)

    #initialization for TR
    tRScore = np.zeros((len(query)+1,len_tR+1),dtype=np.int_)
    for j in xrange(1,len_tR+1):
        tRScore[0][j] = -10000000

    tRMaxScore = np.zeros((len(query)+1,len_tR+1),dtype=np.int_)-1
    tR_rowMax = np.zeros((len(query)+1,2),dtype=np.int_)
    dpF.sw_tR_simple(query,tRArray,tRScore,tRMaxScore,len_prefix,preScore,tR_rowMax) #Calling sw_tR_simple function (implemented in Cython)

    #initialization for suffix
    index_i_tR = 0
    Suffix_BestScore = 0
    sufScore = np.zeros((len(query)+1,len_suffix+1),dtype=np.int_)
    for j in xrange(1,len_suffix+1):
        sufScore[0][j] = -10000000

    sufMaxScore = np.zeros((len(query)+1,len_suffix+1),dtype=np.int_)-1
    dpF.sw_suffix(query,sArray,sufScore,sufMaxScore,tR_rowMax) #Calling sw_suffix function (implemented in Cython)
    index_i_tR,Suffix_BestScore,alignSufxQ = tb_suffix(query,sArray,sufScore,sufMaxScore)

    number_tR = 0.0
    naiveTR = 0.0
    tR_maxScore = 0
    repeatAligned = []
    number_tR, repeatAligned, tR_maxScore, tR_exitScore, naiveTR = finalTrace(query,tR,tRScore,tRMaxScore,tR_rowMax,index_i_tR) #Traceback from the suffix matrix, S and into TR repeat matrix, R 

    alignPrfxQ = tb_prefix(query,pArray,preScore,preMaxScore,Prefix_BestScoreIndex) #Here goes the prefix matrix traceback. here query is in int, but tR is in string representation
    print "exit finalTrace..."
    return number_tR, repeatAligned, alignSufxQ, alignPrfxQ


if __name__ == '__main__': 
    queryfn = sys.argv[1]
    prefixfn = sys.argv[2]
    suffixfn = sys.argv[3]
    trfn = sys.argv[4]
    
    query = list(SeqIO.parse(queryfn, "fasta"))
    prefix = list(SeqIO.parse(prefixfn, "fasta"))
    suffix = list(SeqIO.parse(suffixfn, "fasta"))
    tr = list(SeqIO.parse(trfn, "fasta"))

    assert len(query) == len(prefix)
    assert len(query) == len(suffix)
    assert len(query) == len(tr)

    for i,(q,p,s,t) in enumerate(zip(query,prefix,suffix,tr)):
        q_seq = map(ord,q.seq)
        q_seq = np.asarray(q_seq)
        p_s_t = [p.seq,s.seq,t.seq]
        number_of_tr, repeat_aligned, suffix_aligned, prefix_aligned = alignRegions(q_seq,p_s_t) 
        out = [i,q.id,t.id,t.seq,number_of_tr,
               repeat_aligned[0],repeat_aligned[1],
               suffix_aligned[0],suffix_aligned[1],
               prefix_aligned[0],prefix_aligned[1]]
        print "\t".join(map(str,out))
