#!/bin/env python
import numpy as np

def smith_waterman_prefix(query,prefix):
    query_size = len(query)
    prefix_size = len(prefix)

    prefix_score_matrix = np.zeros((query_size,prefix_size),dtype=np.int_)
    for j in range(1,prefix_size+1):
        prefix_score_matrix[0][j] = -5*j

    mismatch = -6
    match = 5
    penalty = -5
    possible_vals = [0,0,0]
    for i in range(1,query_size+1):
        for j in range(1,prefix_size+1):            
            if query[i-1] == prefix[j-1]:
                diag = match
            else:
                diag = mismatch
            possible_vals[0] = prefix_score_matrix[i-1,j-1] + diag
            possible_vals[1] = prefix_score_matrix[i-1,j] + penalty
            possible_vals[2] = prefix_score_matrix[i,j-1] + penalty
            max_score = possible_vals[0]
            max_score_index = 0
            for k in [1,2]:
                if possible_vals[k] > max_score:
                    max_score = possible_vals[k]
                    max_score_index = k
            prefix_score_matrix[i,j] = max_score
    
    best_score = 0
    column_index = prefix_size
    for i in range(query_size + 1):
        if best_score < prefix_score_matrix[i,column_index]:
            best_score = prefix_score_matrix[i,column_index]
            column_index = i
    
    return best_score,column_index,prefix_score_matrix

def smith_waterman_tandem_repeat(query,tr,prefix_score_matrix,prefix_size):
    query_size = len(query)
    tr_size = len(tr)

    tr_score_matrix = np.zeros((query_size,tr_size),dtype=np.int_)
    for j in range(1,tr_size+1):
        tr_score_matrix[0][j] = -100000000
        
    mismatch = -6
    match = 5
    penalty = -5
    possible_vals = [0,0,0,0,0]
    for i in range(1,query_size+1):
        for j in range(1,tr_size+1):
            if query[i-1] == tr[j-1]:
                diag = match
            else:
                diag = mismatch
            possible_vals[0] = tr_score_matrix[i-1,j] + penalty
            possible_vals[1] = tr_score_matrix[i,j-1] + penalty
            possible_vals[2] = tr_score_matrix[i-1,j-1] + diag
            possible_vals[3] = tr_score_matrix[i-1,tr_size] + diag + (j-1)*penalty
            possible_vals[4] = prefix_score_matrix[i-1,len_prefix] + diag
            max_score = possible_vals[0]
            max_score_index = 0
            for k in [1,2,3,4]:
                if possible_vals[k] > max_score:
                    max_score = possible_vals[k]
                    max_score_index = k
            tr_score_matrix[i,j] = max_score

def smith_waterman_prefix(query,suffx): 
    query_size = len(query)
    suffix_size = len(suffix)   

    suffix_score_matrix = np.zeros((query_size,suffix_size),dtype=np.int_)
    for j in range(1,suffix_size+1):
        suffix_score_matrix[0][j] = -100000000
    
    mismatch = -6
    match = 5
    penalty = -5
    possible_vals = [0,0,0]
    
    for i in range(1,query_size+1):
        if query[i-1] == suffix[0]:
            diag = match
        else:
            diag = mismatch
        suffix_score_matrix[i-1,1] = 


def run_dp_algorithm(query,prefix,suffix,tr):
    prefix_best_score, prefix_column_index,prefix_score_matrix = smith_waterman_prefix(query,prefix)
    smith_waterman_tandem_repeat(query,tr,prefix_score_matrix,len(prefix))
    smith_waterman_prefix(query,suffx)
