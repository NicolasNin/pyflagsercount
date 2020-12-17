import numpy as np
import pyflagsercount
import pyflagser #to compare result


N=1000
p=0.1
np.random.seed(0) 
m=np.random.random((N,N))<p
np.fill_diagonal(m,False)

m_filtered=np.copy(m)
m_filtered[m!=0]=np.arange(np.count_nonzero(m))
simplex_count=np.array(pyflagser.flagser_count_unweighted(m))
tree_vertex,tree_child=pyflagsercount.flagser_simplex_tree(m)
contain=pyflagsercount.flagser_contain(m)
def all_width(tree_vertex):
    all_len=[]
    for root in  tree_vertex:
        all_len.append(np.array([len(x) for x in root]))
    return pyflagsercount.listoflistToarray(all_len)
all_w=all_width(tree_vertex)
# since all the node of the tree are a simplex we can sum this and get back simplex count
assert (np.all(all_w.sum(axis=0)== simplex_count))
"""
%time pyflagser.flagser_count_unweighted(m)
%time a=pyflagsercount.flagser_contain(m)
%time a=pyflagsercount.flagser_filtered(m_filtered)
%time a=pyflagsercount.flagser_simplex_tree(m)
"""