# Lattice that computes attributes with SAGE

from sage.all import *


class Lattice:

    def __init__(self, gram):
        self.gram = matrix(gram)

    @property
    def center_density(self):
        g = self.gram
        t = g.LLL_gram()
        d = g.nrows()
        reduced_g = t.transpose() * g * t
        r = reduced_g[0][0] ** .5 / 2

        # The following only computes the values necessary to find t.transpose()*g*t[0][0]
        # but this is slower when the matrix size is small:

        # vec = vector([t.column(0) * g.column(i) for i in range(d)])
        # r = (vec * t.column(0)) ** .5 / 2

        return r ** d / g.determinant() ** .5


if __name__ == '__main__':
    L = Lattice([[1, 1, 1], [-1, 0, 2], [3, 5, 6]])
    assert type(L.gram).__name__ == 'Matrix_integer_dense'
    print L.center_density
