array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
       [0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
       [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0, 1, 0, 1, 0, 0],
       [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
       [0, 1, 0, 0, 0, 0, 0, 0, 0, 1],
       [1, 1, 0, 0, 0, 1, 0, 0, 0, 0],
       [0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
       [1, 1, 0, 0, 0, 0, 0, 1, 0, 1]])

[[0 1 0 1 1]
 [0 0 1 1 1]
 [1 1 0 1 1]
 [0 0 1 0 1]
 [1 0 1 1 0]]

 [[5, 15, 24, 16, 2], [0, 0, 0, 7, 2]]
8
0
0
MAXIMAL TRUE
[[1, 2, 3, 4],
[1, 2, 4, 3],
[1, 3, 2, 4],
[1, 3, 4, 2],
[1, 4, 2, 3],
[1, 4, 3, 2],
[2, 4, 0, 3],
[4, 2, 0, 3]])


array([[0, 1, 3, 4],
       [0, 1, 4, 3],
       [1, 2, 3, 4],
       [1, 2, 4, 3],
       [1, 3, 2, 4],
       [1, 3, 4, 2],
       [1, 4, 2, 3],
       [1, 4, 3, 2],
       [2, 0, 1, 3],
       [2, 0, 1, 4],
       [2, 0, 3, 4],
       [2, 0, 4, 3],
       [2, 1, 3, 4],
       [2, 1, 4, 3],
       [2, 4, 0, 3],
       [4, 2, 0, 3]])
array([[2, 0, 1, 3, 4],
       [2, 0, 1, 4, 3]])

import flagcount

import numpy as np
import pyflagser
import pyflagsercount
def unique_identifier(mat):
    dt = np.dtype([('raw', np.void, mat.shape[1] * mat.dtype.itemsize)])
    return mat.reshape(np.prod(mat.shape)).view(dt)


def find_maximal_simplices_from_all(simplices, simplices_higher_dim):
    one_index = unique_identifier(simplices)
    simplices_higher_dim_stacked = unique_identifier(np.vstack([simplices_higher_dim[:, np.delete(np.arange(simplices_higher_dim.shape[1]), x)] for x in range(simplices_higher_dim.shape[1])]))
    return simplices[np.logical_not(np.isin(one_index, simplices_higher_dim_stacked)), :]


all_len=[]
for i in range(10):
     all_len.append(np.array([len(x) for x in a[i]]))
    all_len=pyflagsercount.listoflistToarray(all_len)

m=np.random.random((10,10))<0.5
np.fill_diagonal(m,False)
print(m.astype(int))
print(pyflagser.flagser_count_unweighted(m))
_,ls=flagcount.flagser(m,return_simplices=True)
print(pyflagsercount.flagser_maximal(m))
for i in range (len(ls)-1):
    print(len(find_maximal_simplices_from_all(ls[-2-i],ls[-1-i])))

import scipy.sparse as sparse
sparse.coo_matrix(m).row
sparse.coo_matrix(m).col


[0 0 1 1 1 2 2 2 3]
[1 2 0 2 3 0 1 3 0]

[[4, 9, 10, 2], [0, 0, 6, 2]] should be 0 0 4 2
[[0, 1],
[0, 2],
[1, 0],
[1, 2],
[1, 3],
[2, 0],
[2, 1],
[2, 3],
[3, 0]]), 

[0, 1, 2],  *
[0, 2, 1],  *
[1, 0, 2], *
    [1, 2, 0], *
    [1, 2, 3],
    [1, 3, 0],
[2, 0, 1],
    [2, 1, 0],
    [2, 1, 3],
    [2, 3, 0]

[1, 2, 3, 0],
[2, 1, 3, 0]




array([0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 4, 4], dtype=int32)
array([1, 2, 3, 4, 0, 3, 4, 0, 1, 0, 1, 2], dtype=int32)
[[5, 12, 11, 1], [0, 0, 8, 1]]
                      0 7

[[0, 1],
[0, 2],
[1, 0],
[1, 2],
[1, 3],
[2, 0],
[2, 1],
[2, 3],
[3, 0]
[0, 1, 2],
[0, 2, 1],
[1, 0, 2],
[1, 2, 0],
[1, 2, 3],
[1, 3, 0],
[2, 0, 1],
[2, 1, 0],
[2, 1, 3],
[2, 3, 0]]
[1, 2, 3, 0],
[2, 1, 3, 0]]
