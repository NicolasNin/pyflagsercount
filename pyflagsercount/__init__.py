"""Implementation of the python API for the cell count of the flagser C++ library."""

import numpy as np
import scipy.sparse as sparse
from .pyflagsercount import compute_cell_count,compute_cell_count_filtered,compute_cell_count_maximal,compute_simplex_tree

def listoflistToarray(l):
    """ convert a list of list of number as a numpy array, fill with 0 when a list has not the same size"""
    max_dim=max([len(c) for c in l])
    all_array=[np.pad(c,(0,max_dim-len(c)),"constant",constant_values=(0,0)) for c in l]
    return np.array(all_array)

def convertCOO(adjacency_matrix, ret_data=True):
    """ simple utiles to convert array or sparse to row col value """
    if isinstance(adjacency_matrix, np.ndarray) :
        row,col = np.nonzero(adjacency_matrix)
        value = adjacency_matrix[adjacency_matrix != 0]
    elif isinstance(adjacency_matrix, sparse.spmatrix):
        if not isinstance(adjacency_matrix,sparse.coo_matrix):
            adjacency_matrix=adjacency_matrix.tocoo()
        row = adjacency_matrix.row
        col = adjacency_matrix.col
        value =adjacency_matrix.data

    if ret_data:
        return row,col,value
    else:
        return row,col
def flagser_simplex_tree(adjacency_matrix):
    """
    Return the simplex tree 
    :param adjacency_matrix: either a numpy 2d square-array with *positive* value for filtration, or a scipy.sparse
    matrix (will be converted to sparse coo matrix in any case)
    :return: contains values for each vertex and each dimension
    """
    N=adjacency_matrix.shape[0]
    row,col=convertCOO(adjacency_matrix,ret_data=False)
    return compute_simplex_tree(N, np.transpose(np.array( (row,col))))

def flagser_contain(adjacency_matrix):
    """
    Compute for each vertex, the number of times this vertex appears in a simplex for each dimension
    :param adjacency_matrix: either a numpy 2d square-array with *positive* value for filtration, or a scipy.sparse
    matrix (will be converted to sparse coo matrix in any case)
    :return: contains values for each vertex and each dimension
    """
    N=adjacency_matrix.shape[0]
    row,col=convertCOO(adjacency_matrix,ret_data=False)
    return compute_cell_count(N, np.transpose(np.array( (row,col))))

def flagser_maximal(adjacency_matrix):
    N=adjacency_matrix.shape[0]
    row,col=convertCOO(adjacency_matrix,ret_data=False)
    return compute_cell_count_maximal(N, np.transpose(np.array( (row,col))))

def flagser_filtered(adjacency_matrix, return_as_array=True):
    """
    Compute simplex count over a filtration using max value
    :param adjacency_matrix: either a numpy 2d square-array with *positive INTEGER* value for filtration, or a scipy.sparse
    matrix (still *POSITIVE INTEGER VALUE* )
    :param return_as_array: boolean , if true the output will be a numpy array of shape (filtration_size,maximal dimension) 
    :return: return simplex count over a filtration as a list of list or numpy array
    """

    N=adjacency_matrix.shape[0]
    if isinstance(adjacency_matrix, np.ndarray) :
        row,col = np.nonzero(adjacency_matrix)
        value = adjacency_matrix[adjacency_matrix != 0]
    elif isinstance(adjacency_matrix, sparse.spmatrix):
        if not isinstance(adjacency_matrix,sparse.coo_matrix):
            adjacency_matrix=adjacency_matrix.tocoo()
        row = adjacency_matrix.row
        col = adjacency_matrix.col
        value =adjacency_matrix.data

    if return_as_array:
        cc=compute_cell_count_filtered(adjacency_matrix.shape[0], np.transpose(np.array((row,col,value))))
        return listoflistToarray(cc)
    else:
        return compute_cell_count_filtered(adjacency_matrix.shape[0], np.transpose(np.array((row,col,value))))

