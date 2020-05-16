# from functools import lru_cache


class Vector(list):

    def __init__(self, *args):
        if len(args) == 1:
            super(Vector, self).__init__(*args)
        elif len(args) > 1:
            super(Vector, self).__init__(list(args))

    # @lru_cache(maxsize=None)
    def sdot(self):
        return self.dot(self)

    # @lru_cache(maxsize=None)
    def dot(self, rhs):
        rhs = Vector(rhs)
        assert len(self) == len(rhs)
        return sum(map(lambda x: x[0] * x[1], zip(self, rhs)))

    # @lru_cache(maxsize=None)
    def proj_coeff(self, rhs):
        rhs = Vector(rhs)
        assert len(self) == len(rhs)
        return self.dot(rhs) / self.sdot()

    # @lru_cache(maxsize=None)
    def proj(self, rhs):
        rhs = Vector(rhs)
        assert len(self) == len(rhs)
        return self.proj_coeff(rhs) * self

    # @lru_cache(maxsize=None)
    def __add__(self, other):
        assert len(self) == len(other)
        return Vector(map(sum, zip(self, other)))

    def __sub__(self, rhs):
        rhs = Vector(rhs)
        assert len(self) == len(rhs)
        return Vector(x - y for x, y in zip(self, rhs))

    def __mul__(self, rhs):
        if type(rhs) == Vector:
            return self.dot(rhs)
        return Vector(x * rhs for x in self)

    def __rmul__(self, lhs):
        return Vector(x * lhs for x in self)

    def __repr__(self):
        return "[{}]".format(", ".join(str(x) for x in self))

    def __hash__(self):
        return hash(repr(self))


if __name__ == '__main__':
    pass
