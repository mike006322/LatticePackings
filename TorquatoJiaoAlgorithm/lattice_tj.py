"""
Torquato-Jiao (TJ) algorithm for finding sphere packings with high center density
Adapted for lattice packings in the following paper:
https://arxiv.org/pdf/1304.5003.pdf
"""
import logging
from core.log_util import log
from core.lll import lll_reduction
from core.norms import euclidean_norm as norm, sum_of_squared_coefficietns
from core.lattice_enumeration import find_vectors_less_than
from simplex_method import *
from core.matrix import Matrix
from core.polynomial import Polynomial
import numpy as np
from core.lattice import Lattice


@log
def lattice_tj(m):
    """
    input m, that represents a basis for a lattice sphere packing
    returns lattice packing with a density that is a local maximum
    """
    n = len(m)
    b = np.array(lll_reduction(m, 0.75))  # perform LLL reduction on the matrix
    D = norm(b[0])  # D is length of shortest vector because b is LLL reduced
    R = D * 1.1
    epsilon = Matrix.identity(n)
    threshold = 1e-12
    while sum_of_squared_coefficietns(epsilon) > threshold:
        det_b = np.linalg.det(b)
        shortest_vectors = find_vectors_less_than(b.transpose(), R)
        shortest_vectors = remove_zero_vector(shortest_vectors)
        constraints = make_constraints(shortest_vectors, D, R, n)
        epsilon = make_epsilon(constraints, n)
        updated_b = b + b @ epsilon
        updated_det = np.linalg.det(updated_b)
        # if det_b > 0:
        #     if updated_det > det_b:
        #         logging.debug('Updated determinant large. Halving epsilon')
        #         epsilon *= .5
        #         updated_b = b + epsilon @ b
        #         logging.debug('updated determinant ' + str(updated_det))
        #         updated_det = np.linalg.det(updated_b)
        # else:
        #     if updated_det < det_b:
        #         logging.debug('Updated determinant large. Halving epsilon')
        #         epsilon *= .5
        #         updated_b = b + epsilon @ b
        #         logging.debug('updated determinant ' + str(updated_det))
        #         updated_det = np.linalg.det(updated_b)
        b = updated_b
        print('center density = ' + str(Lattice(b).center_density))

    return b


def make_epsilon(constraints, n):
    simplex_input = make_simplex_input(constraints, n)
    simplex_input = np.array(simplex_input)
    epsilon_variables = log(simplex_method_scipy, show_input=False)(simplex_input, unrestricted=True)
    epsilon = np.array(reshape_epsilon(epsilon_variables, n))
    return epsilon


def remove_zero_vector(vectors):
    """
    Removes the all zero vector from list of shortest vectors
    This step is necessary when working with floating point numbers because often 0 is represented as 1e-x,
    where x is large, and later it will add an unwanted constraint in the optimization
    """
    assert type(vectors[0]).__name__ == 'ndarray'
    res = vectors.copy()
    zero = np.zeros(len(vectors[0]))
    ind = 0
    size = len(vectors)
    while ind != size and not np.array_equal(res[ind], zero):
        ind += 1
    if ind != size:
        res.pop(ind)
    return res


def make_constraints(shortest_vectors, D, R_i, n):
    """
    makes constraints out of the shortest vectors
    as per algorithm specifications
    to be used in simplex method
    """
    epsilon = np.zeros((n, n)).tolist()
    epsilon = Matrix(epsilon)
    variable_numbers = []
    # fill epsilon with variables
    variable_index = 0
    for i in range(len(epsilon)):
        for j in range(len(epsilon[0])):
            if i > j:
                epsilon[i][j] = epsilon[j][i]
            else:
                epsilon[i][j] = Polynomial('x' + str(variable_index))
                variable_index += 1
                variable_numbers.append((i, j))

    constraints = []
    # shortest vector constraints
    if shortest_vectors:
        for vector in shortest_vectors:
            vector = Matrix([vector])
            constraint = -1 * vector * epsilon * vector.transpose()
            constraint = make_vector_from_linear_polynomial(constraint[0][0], n)
            v_v_t = vector * vector.transpose()
            v_v_t = v_v_t[0][0]  # change from Matrix type to just a number
            constraint.append(-1 * (D ** 2 - v_v_t) / 2)
            constraints.append(constraint)
    # example [1, 2, 3, 2, 4, 6, 3, 6, 9, -5]
    # meaning x_0 + 2x_1 + 3x_2 + 2x_3 + 4x_4 + 5x_5 + 3x_6 + 6x_7 + 9x_8 >= -5

    lam = (1 - (D / R_i) ** 2) / 2  # lambda
    # bound lowest eigenvalue of epsilon from below by -lamda
    # -.5*lam <= diagonal element of epsilon
    # -.5*lam/(d-1) <= off-diagonal element of epsilon
    # off-diagonal element of epsilon <= .5*lam(d-1)
    # making this all into form 'element of epsilon <= #':
    #
    # - diagonal element of epsilon <= .5*lam
    # - off-diagonal element of epsilon <= .5*lam/(d-1)
    # off-diagonal element of epsilon <= .5*lam(d-1)

    for v, variable in enumerate(variable_numbers):
        i, j = variable
        if i == j:  # if diagonal element of epsilon
            # - diagonal element of epsilon <= .5*lam
            constraint = []
            for _ in range(len(variable_numbers)):
                constraint.append(0)
            constraint[v] = -1
            constraint.append(.5 * lam)
            constraints.append(constraint)
        else:  # if off-diagonal element of epsilon
            # -off-diagonal element of epsilon <= .5*lam/(d-1)
            # off-diagonal element of epsilon <= .5*lam(d-1)
            constraint_positive = []
            constraint_negative = []
            for _ in range(len(variable_numbers)):
                constraint_positive.append(0)
                constraint_negative.append(0)
            constraint_positive[v] = 1
            constraint_negative[v] = -1
            constraint_positive.append(.5 * lam / (n - 1))
            constraint_negative.append(.5 * lam / (n - 1))
            constraints.append(constraint_positive)
            constraints.append(constraint_negative)

    return constraints


def make_simplex_input(constraints, n):
    """
    we need to multiply the constraints by negative one because
    "simplex_method" is standardized to "<=" vectors instead of ">=" vectors

    output matrix:
    constraints         [constraint vars ][constants ]
    .                   .                   .
    .                   .                   .
    .                   .                   .

    objective function  [obj func vars  ][          ]

    example:
    epsilon =
    [[2, 3]
     [2, 3]]
    objective function: -2 + -3 <= 0  -> [-2, -3, 0]
    constraints:
    4x_0 + 3x_1 <= 7  -> [4, 3, 7]
    -13x_0 <= 2  ->  [-13, 0, 2]
    output matrix:
    [  4, 3, 1, 0, 0, 7]
    [-13, 0, 0, 1, 0, 2]
    [-2, -3, 0, 0, 1, 0]
    """
    # make the objective function vector
    objective_function = []
    for i in range(n):
        for j in range(i, n):
            if i == j:
                objective_function.append(1)
                # objective_function.append(epsilon[i][j])
            else:
                objective_function.append(0)
    objective_function.append(0)

    # begin to make the output matrix by putting the objective function after the constraints

    output_matrix = Matrix(constraints)
    output_matrix.append(objective_function)
    # output_matrix now has length n**2 + 1, where epsilon is n by n
    return output_matrix


def make_vector_from_linear_polynomial(poly, n):
    """
    input poly, Polynomial
    n, length of vector
    2x_1 + 3x_2 + x_3 + 8x_4, 4  ->  [2, 3, 1, 8]
    5x_2 + 3x_4, 5 ->  [0, 0, 0, 3, 0]
    """
    number_of_variables = n * (n + 1) // 2
    # initialize res to be all zero's
    res = []
    for i in range(number_of_variables):
        res.append(0)

    # add non-linear terms to make sure all variables are present
    for i in range(number_of_variables):
        poly += Polynomial('x' + str(i) + '^2')

    for term in poly.term_matrix[1:]:
        for variable in range(1, number_of_variables + 1):
            if term[variable] == 1:
                res[variable - 1] = term[0]
                break
    return res


def reshape_epsilon(variables, n):
    epsilon = []
    variable_index = 0
    for i in range(n):
        epsilon.append([])
        for j in range(n):
            if i > j:
                epsilon[i].append(epsilon[j][i])
            else:
                epsilon[i].append(variables[variable_index])
                variable_index += 1
    return epsilon


def main():
    logging.basicConfig(filename='lattice_tj.log',
                        level=logging.DEBUG,
                        format='%(asctime)s - %(name)s - %(threadName)s -  %(levelname)s - %(message)s',
                        filemode='w')
    m = [[1, 12, 1, 3, 0], [-3, 0, 2, 3, 0], [3, 5, 6, 4, 0], [3, -4, -5, 6, 0], [-8, 4, 7, -3, 7]]
    # m = [[1, 1, 1, 3], [-1, 0, 2, 3], [3, 5, 6, 4], [3, -4, -5, 6]]
    # m = [[1, 1, 1], [-1, 0, 2], [3, 5, 6]]
    # m = Matrix.identity(3)
    # m = [[2, 1, 4], [18, -3, 0], [-3, 1, 6]]

    print(Lattice(m).center_density)
    denser_matrix = lattice_tj(m)
    print('Output: \n', denser_matrix)
    print(Lattice(denser_matrix.tolist()).center_density)


def get_density_test():
    m = [[-2992.9880115283113, -10.72655918418694, -373.644364082849, -79.76874811239823],
         [1118.502728914449, 5.872511951906576, 139.03821745089388, 31.520105005997586],
         [1183.3611207871302, 2.8580632518764837, 148.10165826317575, 30.42962084140757],
         [2217.808564341084, 8.538061925470647, 276.7082101552752, 59.66740931905222]]
    # actual center density: 0.0593872379809836881793289678657
    # My algorith, said 0.13181955941418586
    # (4, [-2992.9880115283113, -10.72655918418694, -373.644364082849, -79.76874811239823, 1118.502728914449, 5.872511951906576, 139.03821745089388, 31.520105005997586, 1183.3611207871302, 2.8580632518764837, 148.10165826317575, 30.42962084140757, 2217.808564341084, 8.538061925470647, 276.7082101552752, 59.66740931905222]])
    print(Lattice(m).center_density)


if __name__ == '__main__':
    main()
    # get_density_test()
