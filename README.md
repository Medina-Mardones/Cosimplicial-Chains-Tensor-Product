# Constructing a counterexample for the CPC-conjecture 


Consider the tensor product of the chains of a d-simplex with themselves. 

Bdom be the canonical ordered basis for the degree d-part of this complex

print('\nConsider the involution T given by transposition of factors. It induces an involution on B.')
print('The following list contains the positions, with respect to the order of B, of transpositionally related elements (b, Tb):\n')
pairs = get_transp_pairs(dim)
print(pairs)
print('\nA vector $c=(c_i)$ with respect to this basis is called unhampered if for each $i$ the product $c_i (Tc)_i$ is 0.')

print('\nLet Bcod be the following ordered basis for the degree',dim-1,'part of this complex:\n')
Bcod = get_tensor_gr_gen(dim)[dim-1]
print(Bcod, len(Bcod))

print('\nNow consider the kernel K of the matrix M representing the boundary in degree',dim,'with respect to B and B\' \n')
print(get_tensor_partial(dim))
        
print('\nAn element $c$ is called cpc if $(s_i \otimes s_i)(c)$ is 0 for each $i$. Here $s_i$ is the induced $i$-degeracy map.\n')
print('Conjecture: The only unhampered vector representing a cpc element in K is 0.\n')
