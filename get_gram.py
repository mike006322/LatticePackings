# translates from a given search matrix to a gram matrix

from search_utilities import *
from sage_lattice import *
from matrix import Matrix
from trace_form_symmetric_matrices import *


def search_to_gram(p, f, search_matrix, g, m):
    """
    input: p - prime dimension,
    f - conductor,
    search_matrix,
    g - symmetric matrix of trace form,
    m - modulo factor
    """
    zeta = find_primitive_root(p, f)
    h = make_h(zeta, p, f)
    h = m * h
    j = f * Matrix(search_matrix)
    t = h.concatenate(j)
    t = augment_identity_times_factor(t, m * f)
    t = matrix(t)  # make t a SAGE matrix
    n = get_nullspace(t)
    n = n.matrix_from_columns(range(p))  # SAGE specific. Python alternative: n = Matrix(n[:p]).transpose()
    gram_matrix = n * g * n.transpose()
    d = get_density(gram_matrix)
    return gram_matrix, d


def get_density(gram_matrix):
    l = Lattice(gram_matrix)
    return l.center_density


def get_nullspace(m):
    # if using Matrix class from matrix.py: return m.nullspace()
    return m.transpose().kernel().basis_matrix()  # SAGE specific


def dim_13():
    p = 13
    f = 53
    g = matrix(DIM_13_TR_SYM_MATRIX)
    m = 2

    # top_part_right = search_matrix
    # id_4 = Matrix.identity(4)
    # top_part_left = id_4.concatenate(Matrix([[0], [0], [0], [0]]), axis=1)
    # top_part = top_part_left.concatenate(top_part_right, axis=1)
    # bottom_part = Matrix([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])
    # search_matrix = top_part.concatenate(bottom_part)

    search_matrix = [[1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1],
                     [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
                     [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0],
                     [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1],
                     [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]

    search_matrix = [[1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1],
                     [0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0],
                     [0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0],
                     [0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1],
                     [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]

    return search_to_gram(p, f, search_matrix, g, m)


def dim_11(search_matrix):
    p = 11
    f = 23
    m = 2
    g = matrix(DIM_11_TR_SYM_MATRIX)
    zeta = find_primitive_root(p, f)
    h = m * make_h(zeta, p)
    top_part_right = search_matrix
    id_4 = Matrix.identity(4)
    top_part_left = id_4.concatenate(Matrix([[0], [0], [0], [0]]), axis=1)
    top_part = top_part_left.concatenate(top_part_right, axis=1)
    bottom_part = Matrix([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])
    search_matrix = top_part.concatenate(bottom_part)
    j = f * search_matrix
    t = h.concatenate(j)
    t = augment_identity_times_factor(t, 2 * f)
    t = matrix(t)  # make t a SAGE matrix
    n = get_nullspace(t)
    n = n.matrix_from_columns(range(p))  # SAGE specific. Python alternative: n = Matrix(n[:p]).transpose()
    gram_matrix = n * g * n.transpose()
    print(gram_matrix.rank())
    print(gram_matrix.transpose().rank())
    print(g.is_positive_definite())
    d = get_density(gram_matrix)

    return search_to_gram(p, f, search_matrix, g, m)


def dim_11_high(search_matrix):
    p = 11
    f = 23
    m = 2
    g = matrix(DIM_11_TR_SYM_MATRIX)
    zeta = find_primitive_root(p, f)
    h = m * make_h(zeta, p)
    top_part_right = search_matrix
    id_9 = Matrix.identity(9)
    top_part_left = id_9.column_sub_matrix(8)
    # top_part_left = id_9.concatenate(Matrix([[0], [0], [0], [0], [0], [0], [0], [0], [0]]), axis=1)
    top_part = top_part_left.concatenate(top_part_right, axis=1)
    bottom_part = Matrix([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])
    search_matrix = top_part.concatenate(bottom_part)
    j = f * search_matrix
    t = h.concatenate(j)
    t = augment_identity_times_factor(t, 2 * f)
    t = matrix(t)  # make t a SAGE matrix
    n = get_nullspace(t)
    n = n.matrix_from_columns(range(p))  # SAGE specific. Python alternative: n = Matrix(n[:p]).transpose()
    gram_matrix = n * g * n.transpose()
    d = get_density(gram_matrix)

    return search_to_gram(p, f, search_matrix, g, m)


def dim_7(search_matrix):
    p = 7
    f = 29
    m = 2
    zeta = find_primitive_root(p, f)
    g = matrix(DIM_7_TR_SYM_MATRIX)
    h = m * make_h(zeta, p)
    top_part = search_matrix
    bottom_part = Matrix([[1, 1, 1, 1, 1, 1, 1]])
    search_matrix = top_part.concatenate(bottom_part)
    j = f * search_matrix
    t = h.concatenate(j)
    t = augment_identity_times_factor(t, 2 * f)
    t = matrix(t)  # make t a SAGE matrix
    n = get_nullspace(t)
    n = n.matrix_from_columns(range(p))  # SAGE specific. Python alternative: n = Matrix(n[:p]).transpose()
    gram_matrix = n * g * n.transpose()
    d = get_density(gram_matrix)

    return search_to_gram(p, f, search_matrix, g, m)


def dim_5():
    search_matrix = Matrix([[1, 1, 1, 1, 1]])
    p = 5
    f = 11
    m = 2
    g = matrix(DIM_5_TR_SYM_MATRIX)
    zeta = find_primitive_root(p, f)
    h = m * make_h(zeta, p)
    j = f * search_matrix
    t = h.concatenate(j)
    t = augment_identity_times_factor(t, 2 * f)
    t = matrix(t)  # make t a SAGE matrix
    n = get_nullspace(t)
    n = n.matrix_from_columns(range(p))  # SAGE specific. Python alternative: n = Matrix(n[:p]).transpose()
    gram_matrix = n * g * n.transpose()
    d = get_density(gram_matrix)

    return search_to_gram(p, f, search_matrix, g, m)


def dim_3():
    search_matrix = Matrix([[1, 1, 1]])
    p = 3
    f = 7
    m = 2
    g = matrix(DIM_3_TR_SYM_MATRIX)
    zeta = find_primitive_root(p, f)
    h = m * make_h(zeta, p)
    j = f * search_matrix
    t = h.concatenate(j)
    t = augment_identity_times_factor(t, 2 * f)
    t = matrix(t)  # make t a SAGE matrix
    n = get_nullspace(t)
    n = n.matrix_from_columns(range(p))  # SAGE specific. Python alternative: n = Matrix(n[:p]).transpose()
    gram_matrix = n * g * n.transpose()
    d = get_density(gram_matrix)

    return search_to_gram(p, f, search_matrix, g, m)


if __name__ == '__main__':
    # search_matrix13 = int_to_bi_matrix(402666846, 4, 8)
    # print(dim_13(search_matrix13)[0])

    # search_matrix11 = int_to_bi_matrix(2622716, 4, 6)
    # search_matrix11 = int_to_bi_matrix(16763378, 4, 6)
    # print(search_matrix11)
    # print(list(dim_11(search_matrix11)[0]))
    # print(dim_11(search_matrix11))

    # search_matrix7 = int_to_bi_matrix(410538, 3, 7)
    # print(search_matrix7)
    # print(dim_7(search_matrix7))

    # print(dim_5())

    # search_matrix_11high = int_to_bi_matrix(12900842, 9, 3)
    # print(search_matrix_11high)
    # print(dim_11_high(search_matrix_11high))

    print(dim_3())

    pass
"""
Q := Matrix(11, [ 67712,  29095,  37030,  49197,  13754, -12696,  31740,   9522,  -6348,  30682,  14812, 29095,  37030,  39146,  50255,  37030,   7406,   7406,   8464,  15870,   4232,  23276, 37030,  39146,  48668,  55016,  31740,   2116,  21160,   8464,  10580,  -4232,  29624, 49197,  50255,  55016,  87814,  56074,   9522,  24334, -13754,  26450,  14812,  30682, 13754,  37030,  31740,  56074,  52900,  21160,  -2116,  -4232,  31740,   6348,  14812, -12696,   7406,   2116,   9522,  21160,  44436,  -4232,  -4232,  -4232,  -4232,  -4232, 31740,  7406,  21160,  24334,  -2116,  -4232,  44436,  -4232,  -4232,  -4232,  -4232,  9522,   8464,  8464,  -13754, -4232,  -4232,  -4232, 44436,  -4232,  -4232,  -4232, -6348,  15870,  10580,  26450,  31740,  -4232,  -4232,  -4232,  44436,  -4232,  -4232, 30682,   4232,  -4232,  14812,   6348,  -4232,  -4232,  -4232,  -4232,  44436,  -4232, 14812,  23276,  29624,  30682,  14812,  -4232,  -4232,  -4232,  -4232,  -4232,  44436]);

"""
