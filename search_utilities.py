# function used to search

from matrix import Matrix
from time import time
import logging


def make_h(zeta, p, f=0):
    """
    p, prime dimension
    zeta, primitive p-th root of unity
    n, conductor
    matrix H_1
    """
    # make (p-1)/2 x p  matrix of all 1's
    # don't touch the first column and first row of 1's
    # of the remaining matrix,
    # i = row number, j = column number
    # multiply each entry by zeta^(i-1)(j-1)
    row_size = (p - 1) // 2 + 1
    column_size = p
    h = Matrix.ones((row_size, column_size))
    # h_1 = np.ones([row_size, column_size])
    for i in range(1, row_size):
        for j in range(1, column_size):
            h[i][j] *= zeta ** (i * j)
            if f != 0:
                h[i][j] %= f
    return h


def find_primitive_root(p, n):
    """
    find a pth-primtive root of unity mod n
    """
    for i in range(2, n):
        order = 1
        multiplied = i
        while multiplied % n != 1:
            multiplied *= i
            order += 1
        if order == p:
            return i


def augment_identity_times_factor(matrix, factor):
    """
    augments factor*identity right of matrix
    """
    id = Matrix.identity(len(matrix))
    return Matrix.concatenate(matrix, id*factor, axis=1)


def int_to_bi_matrix(number, n, m):
    """
    Input number, an integer; and input the dimensions of matrix: n, m
    Returns an n by m binary matrix whose entries are the binary form of input number
    """
    res = []
    for i in range(n):
        res.append([])
        for j in range(m):
            if len(bin(number)) - 2 > m * i + j:
                res[i].append(int(bin(number)[len(bin(number)) - 1 - (m * i + j)]))
            else:
                res[i].append(0)
    return Matrix(res)


def add_to_file(d, filename):
    data_file = open(filename, 'a')
    data_file.write(d + '\n')


def configure_log(filename):
    log_filename = 'logs/' + filename + '.log'
    logging.basicConfig(filename=log_filename,
                        level=logging.DEBUG,
                        format='%(asctime)s - %(name)s - %(threadName)s -  %(levelname)s - %(message)s',
                        filemode='w')


def log_progress(current_number, search_size, start_time):
    """
    every 25000 make a long entry
    """
    if current_number % 25000 == 0:
        current_time = time()
        percent_complete = float(current_number) / float(search_size) * 100
        elapsed_time = current_time - start_time
        logging.info(
            'current number {}, percent complete: {:.2f}%, time elapsed: {:.0f} seconds'.format(current_number,
                                                                                                  percent_complete,
                                                                                                  elapsed_time))


def count_lines(filename):
    """
    returns the number of lines in the file
    excludes last line if it's blank
    """
    data_file = open(filename, 'r')
    data_list = data_file.readlines()
    number_of_lines = len(data_list)
    if data_list[-1] == '\n':
        number_of_lines -= 1
    return number_of_lines


if __name__ == '__main__':
    test_matrix_1 = Matrix.identity(2)
    assert list(augment_identity_times_factor(test_matrix_1, 2)) == [[1, 0, 2, 0], [0, 1, 0, 2]]

    # print(int_to_bi_matrix(402666846, 4, 8))

    print(find_primitive_root(13, 53))

    # print(count_lines('search_data/search_dim_11save.txt'))
