"""
Searches dimension 5 for lattice with density 1/(8sqrt(2)) = 0.08838834764
"""

from search_utilities import *
from sage_lattice import *
from matrix import Matrix
from trace_form_symmetric_matrices import DIM_3_TR_SYM_MATRIX


def make_search_matrix(i):
    """
    returns matrix [X, X, X] where X's are the binary digits of integer i
    The index is 2^1, therefore we have only one search row with p columns.
    """
    search_matrix = int_to_bi_matrix(i, 1, 3)
    # search_matrix = search_matrix.concatenate(matrix([1, 1, 1]))
    return Matrix(search_matrix)


def get_density(gram_matrix):
    l = Lattice(gram_matrix)
    return l.center_density


def get_nullspace(m):
    # if using Matrix class from matrix.py: return m.nullspace()
    return m.transpose().kernel().basis_matrix()  # SAGE specific


def main():
    name = 'search_dim_3'
    data_filename = 'search_data/' + name + '.txt'
    p = 3  # prime dimension
    f = 7  # conductor
    zeta = find_primitive_root(p, f)
    h = 2 * make_h(zeta, p)
    # g is symmetric matrix of trace form
    g = matrix(DIM_3_TR_SYM_MATRIX)
    for i in range(2 ** 3):
        j = f * make_search_matrix(i)
        t = h.concatenate(j)
        t = augment_identity_times_factor(t, 2 * f)
        t = matrix(t)  # make t a SAGE matrix
        n = get_nullspace(t)
        n = n.matrix_from_columns(range(p))  # SAGE specific. Python alternative: n = Matrix(n[:p]).transpose()
        gram_matrix = n * g * n.transpose()
        d = get_density(gram_matrix)
        add_to_file(str(d), data_filename)


if __name__ == '__main__':
    main()
