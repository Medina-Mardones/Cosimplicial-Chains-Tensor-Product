#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 11:45:12 2018

@author: amedinam
"""

import numpy as np
import itertools


# Constructing the ordered bases 

def get_gr_gen(dim):
    '''returns the list whose d-entry is the lexicographically ordered list of 
    d-faces of the (dim)-simplex'''  

    gr_gen = []
    for i in range(dim+1):
        els = [list(x) for x in itertools.combinations(set(range(dim+1)), i+1)]
        gr_gen.append(els)
    
    return gr_gen

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

def get_transp_pairs(dim):
    '''returns the set containing the list with pairs of generators that are equal up to 
    transposition'''
    
    gens = get_tensor_gr_gen(dim)[dim]
    transp_pairs = []
    for gen in gens:
        neg=[gen[1],gen[0]]
        for gen2 in gens:
            if gen2 == neg and gens.index(gen) <= gens.index(gen2):
                transp_pairs.append((gens.index(gen),gens.index(gen2)))    

    return transp_pairs


# Constructing boundary matrices    

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


# Stating the conjecture 

dim = 1 

print('\nConsider the tensor product of the chains of a '+str(dim)+'-simplex with themselves.')

print('\nLet Bdom be the following ordered basis for the degree '+str(dim)+' part of this complex:\n')
Bdom = get_tensor_gr_gen(dim)[dim]
print(Bdom )

print('\nConsider the involution T given by transposition of factors. It induces an involution on B.')
print('The following list contains the positions, with respect to the order of B, of transpositionally related elements (b, Tb):\n')
print(get_transp_pairs(dim))
print('\nA vector $c=(c_i)$ with respect to this basis is called unhampered if for each $i$ the product $c_i (Tc)_i$ is 0.')

print('\nLet Bcod be the following ordered basis for the degree',dim-1,'part of this complex:\n')
Bcod = get_tensor_gr_gen(dim)[dim-1]
print(Bcod, len(Bcod))

print('\nNow consider the kernel K of the matrix M representing the boundary in degree',dim,'with respect to B and B\' \n')
print(get_tensor_partial(dim))
        
print('\nAn element $c$ is called cpc if $(s_i \otimes s_i)(c)$ is 0 for each $i$. Here $s_i$ is the induced $i$-degeracy map.\n')
print('Conjecture: The only unhampered vector representing a cpc element in K is 0.\n')


from sympy import Matrix

M = Matrix(get_tensor_partial(dim))
zeros = Matrix(np.zeros((len(Bcod),1)))

print('The subspace K consists of the elements\n')

sol = M.gauss_jordan_solve(zeros, freevar=True)[0]
free = M.gauss_jordan_solve(zeros, freevar=True)[1]

print(np.array(sol))

print('\nwhere\n')
print(np.array(free))
print('\nare/is the free variable(s)')






    





        
    
