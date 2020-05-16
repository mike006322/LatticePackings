"""
Searches dimension 7 for lattice with density 1/16 = 0.06250
data output file: data/search_dim_7.txt
log file: logs/search_dim_7.log
"""

from search_utilities import *
from sage_lattice import *
from matrix import Matrix
from time import time
from multiprocessing import Pool
from trace_form_symmetric_matrices import DIM_7_TR_SYM_MATRIX


def make_search_matrix(i):
    """
    returns matrix
    [X, X, X, X, X, X, X]
    [X, X, X, X, X, X, X]
    [X, X, X, X, X, X, X]
    [1, 1, 1, 1, 1, 1, 1]

    where X's are the binary digits of integer i
    """
    top_part = int_to_bi_matrix(i, 3, 7)
    bottom_part = Matrix([[1, 1, 1, 1, 1, 1, 1]])
    result = top_part.concatenate(bottom_part)
    return result


def get_density(gram_matrix):
    l = Lattice(gram_matrix)
    return l.center_density


def get_nullspace(m):
    # if using Matrix class from matrix.py: return m.nullspace()
    return m.transpose().kernel().basis_matrix()  # SAGE specific


def subprocess(i):
    """
    subprocess to be run in parallel
    """
    log_progress(i, search_size, start_time)
    j = f * make_search_matrix(i)
    t = h.concatenate(j)
    t = augment_identity_times_factor(t, 2 * f)
    t = matrix(t)  # make t a SAGE matrix
    n = get_nullspace(t)
    n = n.matrix_from_columns(range(p))  # SAGE specific. Python alternative: n = Matrix(n[:p]).transpose()
    gram_matrix = n * g * n.transpose()
    d = get_density(gram_matrix)
    if d == target_density:
        add_to_file(str(d) + ', ' + str(i), data_filename)


p = 7  # prime dimension
name = 'search_dim_' + str(p)
data_filename = 'search_data/' + name + '.txt'
configure_log(name)
target_density = 0.06250
f = 29  # conductor
zeta = find_primitive_root(p, f)
h = 2 * make_h(zeta, p)
# g is symmetric matrix of trace form
g = matrix(DIM_7_TR_SYM_MATRIX)
highest_density = 0  # highest density seen so far. If higher, then record the density
search_size = 2 ** 21  # number of matrices to search
start_time = time()


def main():
    pools = Pool(4)
    results = pools.map(subprocess, range(search_size))
    pools.close()
    pools.join()


if __name__ == '__main__':
    main()
