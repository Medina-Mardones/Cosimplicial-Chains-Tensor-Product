#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 15:31:26 2018

@author: amedinam
"""

import numpy as np
import itertools
from sympy import Matrix

dim = 2

def get_gr_gen(dim):
    '''returns the list whose d-entry is the lexicographically ordered list of 
    d-faces of the (dim)-simplex'''  

    gr_gen = []
    for i in range(dim+1):
        els = [list(x) for x in itertools.combinations(set(range(dim+1)), i+1)]
        gr_gen.append(els)
    
    return gr_gen

print('\ngraded generators for dim =',dim)
for d in range(dim+1):
    print('in degree',d,'\n',get_gr_gen(dim)[d])
    

def get_gr_partial(dim): 
    '''returns the list whose d-entry is the np.array (matrix) representing the 
    boundary map from d to (d-1)-chains of the standard (dim)-simplex with bases 
    lexicographically ordered'''

    gr_gen = get_gr_gen(dim)         
    gr_partial = list()
    for degree in range(dim+1):
        
        col_len = len(gr_gen[degree])
        row_len = len(gr_gen[degree-1])
        partial_sub_degree = np.zeros((row_len,col_len), dtype=int) #initializing matrix
        
        for i in range(col_len): #for all generators in the given degree
            for vertex in range(len(gr_gen[degree][i])): #for all the vertices of a generator
                els = list(gr_gen[degree][i])
                del els[vertex]
                for j in range(row_len): #looking among all generators of lower degree
                    if els == gr_gen[degree-1][j]:
                        partial_sub_degree[j,i]=(-1)**int(vertex)
                        break
        
        gr_partial.append(partial_sub_degree)
    
    return(gr_partial)
    
print('\ngraded partial for dim =',dim)
for d in range(dim+1):
    print('and degree',d,'\n',get_gr_partial(dim)[d])


def get_tensor_gr_gen(dim):        
    '''returns the list whose d-entry is the lexicographically ordered list of pairs of faces 
    of the (dim)-simplex whose dimensions add to d'''      

    gr_gen = get_gr_gen(dim)
    tensor_gr_gen = []    
    for i in range (2*dim+1): #initializing the array
        tensor_gr_gen.append([])            
    for i in range(dim+1): #filling the array
        for j in range(dim+1):
            els = [list(x) for x in itertools.product(gr_gen[i],gr_gen[j])]
            tensor_gr_gen[i+j] += els    
    return tensor_gr_gen

print('\ngraded generators for dim',dim, 'of the tensor product')
for d in range(2*dim+1):
    print('in degree',d,'\n',get_tensor_gr_gen(dim)[d])
    

def get_tensor_partial(dim):
    '''computes the matrix representing the dim-boundary in the tensor product 
    (with respect to the canonical basis lexicografically ordered)'''

    gr_gen = get_gr_gen(dim)
    gr_dim = []
    for l in gr_gen:
        gr_dim.append(len(l))
    
    gr_partial = get_gr_partial(dim)
    
    #constructs the first column
    hor_stack = np.array([0]) 
    for i in range(dim):
        ver_len = gr_dim[i]*gr_dim[dim-1-i]
        block = np.full((ver_len,1),ver_len,dtype=int)
        hor_stack = np.vstack((hor_stack,block))
         
    #constructs and stacks the blocks
    for i in range(dim+1):        
        hor_len = gr_dim[i]*gr_dim[dim-i]
        ver_stack = np.full((1,hor_len),hor_len,dtype=int) #constructs the first row   
        ith_sign = 1
        
        for j in range(dim):
            ver_len = gr_dim[j]*gr_dim[dim-1-j]
    
            if j+1 == i: #dx1:
                identity = np.eye(gr_dim[dim-i],dtype=int)
                block = np.kron(gr_partial[i],identity)
                
            elif dim-j == dim-i: #1xd
                identity = np.eye(gr_dim[i],dtype=int)
                block = np.kron(identity,ith_sign*gr_partial[dim-i])
                
            else: 
                block = np.zeros((ver_len,hor_len),dtype=int)
            
            ver_stack = np.vstack((ver_stack,block))
            ith_sign *= -1
            
        hor_stack = np.hstack((hor_stack,ver_stack))
    
    #deleting first row and column
    hor_stack = np.delete(hor_stack, (0), axis=0)    
    out_matrix = np.delete(hor_stack, (0), axis=1)    

    return out_matrix   

print('\npartial in degree',dim,'\n',get_tensor_partial(dim))
        
def row_reduction(M):
    M = Matrix(M)
    red_M = M.rref()[0]
    pivots = M.rref()[1]
    return(red_M, pivots) 
    
print('\nrow reduction of the deg',dim, 'boundary of the tensor and pivot positions\n', row_reduction(get_tensor_partial(dim)))
    
    
    