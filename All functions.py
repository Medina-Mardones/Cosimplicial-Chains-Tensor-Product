''' All functions developed for the project. Some were not needed 
to find the counterexample
'''

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

print('\ngenerators for dim',dim)
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
        

def get_diag_of_boundary(dim):
    '''returns the vector, with entries: -1, 0 or 1, representing the linear 
    combination obtained by taking boundary of the top dimensional generator 
    and then applying the standard diagonal'''
    
    
    def get_boundary(face):
        '''given a face it returns the lexicographically 
        order list of signed generators appearing in its boundary'''
    
        boundary = []
        ith_sign = 1
        for i in face:
            di_face = list(face)
            di_face.pop(i)
            boundary.append([ith_sign,di_face])
            ith_sign *= (-1)
        
        return(boundary)
   
    def get_diag(sign,face):
        '''given a signed face, [plus/minus 1, [u,v,...,w]], it returns the 
        list of signed generators, [plus/minus 1, [[u,v],[v,...,w]]], appearing 
        in its standard diagonal'''
        
        diag = []
        for i in face:
            diag.append([sign,face[:face.index(i)+1],face[face.index(i):]])
    
        return(diag)      
       
    simplex = list(range(dim+1))    
    boundary = get_boundary(simplex)
    
    diag_partial = []
    for face in boundary: #computes all signed terms in the diagonal of the boundary
        aux = get_diag(face[0],face[1])
        diag_partial.extend(aux)
        
    tensor_gen = get_tensor_gr_gen(dim)[dim-1]    
    out_vector = np.zeros((len(tensor_gen),1),dtype=int)    
    
    for gen in tensor_gen: #constructs the vector representing the linear combination
        for term in diag_partial:
            if gen == term[1:]:
                out_vector[tensor_gen.index(gen)] = term[0]
    
    return out_vector

print('\ndiagonal of the boundary of the top generator\n',get_diag_of_boundary(dim))


def get_relations(dim):
    '''returns a list, each entry corresponds to a column of the non-augmented reduced 
    matrix and it contains a list with the coefficients in terms of free variables and 
    the vector b'''
    
    #getting the augmented matrix
    M = Matrix(get_tensor_partial(dim))
    m = M.shape[1]
    aug_M = M.col_insert(m, Matrix(get_diag_of_boundary(dim)))
    
    #getting the row reduced matrix and the list of pivots
    red_M = aug_M.rref()[0]
    pivots = aug_M.rref()[1]
    
    #a dictionary mapping from the column number of a free variables to their position in the list of free variables
    dict_free_vars = {}
    j = 0
    for i in range(m+1):
        if i not in pivots:
            dict_free_vars[j] = i
            j += 1
    
    #writting the dependance of every column to the free variables with the correct signs
    relations = []
    row = 0
    p = len(dict_free_vars)
    
    for i in range(m): #not including last column
        aux=[]
        if i in pivots:
            for k in range(p-1):
                aux.append((-1)*red_M[i-row,dict_free_vars[k]])
            aux.append(red_M[i-row,dict_free_vars[p-1]])
            
        elif i not in pivots:
            for k in range(p):
                aux.append(0)
            aux[row] = 1
            row += 1
        relations.append(aux)
        
    return relations

print('\nlist of relations\n',get_relations(dim))


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

print('\ncollection of basis elements related by transposition\n',get_transp_pairs(dim))

