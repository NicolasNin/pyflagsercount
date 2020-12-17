import numpy as np
import pyflagsercount
import pyflagser #to compare result


N=2000
p=0.1
m=np.random.random((N,N))<p
np.fill_diagonal(m,False)
simplex_count=np.array(pyflagser.flagser_count_unweighted(m))
tree=pyflagsercount.flagser_simplex_tree(m)
contain=pyflagsercount.flagser_contain(m)
def all_width(tree_vertex):
    all_len=[]
    for root in  tree_vertex:
        all_len.append(np.array([len(x) for x in root]))
    return pyflagsercount.listoflistToarray(all_len)
all_w=all_width(tree)
# since all the node of the tree are a simplex we can sum this and get back simplex count
assert (np.all(all_w.sum(axis=0)== simplex_count))
"""
%time pyflagser.flagser_count_unweighted(m)
%time a=pyflagsercount.flagser_contain(m)
%time a=pyflagsercount.flagser_simplex_tree(m)
"""