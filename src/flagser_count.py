"""Implementation of the python API for the cell count of the flagser C++ library."""

import numpy as np
from pyflagsercontain import compute_cell_count,compute_cell_count_filtered

def flagser_contain(adjacency_matrix):
    return compute_cell_count(adjacency_matrix.shape[0], np.transpose(np.array(np.nonzero(adjacency_matrix))))

def flagser_filtered(adjacency_matrix):
    N=adjacency_matrix.shape[0]
    row,col=np.nonzero(m)
    value=m[m!=0]
    return compute_cell_count_filtered(adjacency_matrix.shape[0], np.transpose(np.array((row,col,value))))

