"""
Searches dimension 11 for lattice with density 1/(16*sqrt(2)) = 0.04419417382
data output file: data/search_dim_11.txt
log file: logs/search_dim_11.log
"""

from search_utilities import *
from sage_lattice import *
from matrix import Matrix
from time import time
from multiprocessing import Pool
from trace_form_symmetric_matrices import DIM_11_TR_SYM_MATRIX


def make_search_matrix(i):
    """
    returns matrix
    [1, 0, 0, 0, 0, 0, 0, 0, X, X, X]
    [0, 1, 0, 0, 0, 0, 0, 0, X, X, X]
    [0, 0, 1, 0, 0, 0, 0, 0, X, X, X]
    [0, 0, 0, 1, 0, 0, 0, 0, X, X, X]
    [0, 0, 0, 0, 1, 0, 0, 0, X, X, X]
    [0, 0, 0, 0, 0, 1, 0, 0, X, X, X]
    [0, 0, 0, 0, 0, 0, 1, 0, X, X, X]
    [0, 0, 0, 0, 0, 0, 0, 1, X, X, X]
    [0, 0, 0, 0, 0, 0, 0, 0, X, X, X]
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    where X's are the binary digits of integer i
    """
    top_part_right = int_to_bi_matrix(i, 9, 3)
    id_9 = Matrix.identity(9)
    top_part_left = id_9.column_sub_matrix(8)
    # top_part_left = id_9.concatenate(Matrix([[0], [0], [0], [0], [0], [0], [0], [0], [0]]), axis=1)
    top_part = top_part_left.concatenate(top_part_right, axis=1)
    bottom_part = Matrix([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])
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
    global max_density_found
    log_progress(i, search_size, start_time)
    j = f * make_search_matrix(i)
    t = h.concatenate(j)
    t = augment_identity_times_factor(t, 2 * f)
    t = matrix(t)  # make t a SAGE matrix
    n = get_nullspace(t)
    n = n.matrix_from_columns(range(p))  # SAGE specific. Python alternative: n = Matrix(n[:p]).transpose()
    gram_matrix = n * g * n.transpose()
    d = get_density(gram_matrix)
    if d > max_density_found:
        max_density_found = d
        add_to_file(str(d) + ', ' + str(i), data_filename)
    if d > target_density:
        add_to_file(str(d) + ', ' + str(i), data_filename)


p = 11  # prime dimension
name = 'search_dim_' + str(p) + '_higher'
data_filename = 'search_data/' + name + '.txt'
configure_log(name)
max_density_found = 0  # to be changed during search
target_density = 0.044
f = 23  # conductor
zeta = find_primitive_root(p, f)
h = 2 * make_h(zeta, p)
# g is symmetric matrix of trace form
g = matrix(DIM_11_TR_SYM_MATRIX)
highest_density = 0  # highest density seen so far. If higher, then record the density
search_size = 2 ** 27  # number of matrices to search
start_time = time()


def main():
    pools = Pool(8)
    results = pools.map(subprocess, range(search_size))
    pools.close()
    pools.join()


if __name__ == '__main__':
    main()
