# testing SAGE on Pycharm python IDE
# Setup instructions are found here:
# https://ask.sagemath.org/question/39742/make-pycharm-recognise-the-sage-python-interpreter/

from sage.all import *

from matrix import *
m = Matrix([[1, 2]])
print m
m = matrix(QQ, 3, 3, [[0, .5, 0], [1, 0, 1], [-1, 0, 2]])
print m.LLL()