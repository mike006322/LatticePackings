# Symmetric matrices of integral trace form of odd prime degree
#
# E. L. d. Oliveira, J. C. Interlando, T. P. d. N. Neto, and J. O. D.Lopes,
# The integral trace form of cyclic extensions of odd prime degree,
# RockyMountain Journal of Mathematics, 47 (2017), pp. 1075-1088.
#
# p = odd prime dimension, f = field conductor, a in Z^n
# Trace form is the quadratic form: f*sum(a_i^2) - ((f-1)/p)sum(a_i)^2
# Below are the symmetric matrices of the quadratic form for different dimensions

DIM_3_TR_SYM_MATRIX = [[5, -2, -2],
                       [-2, 5, -2],
                       [-2, -2, 5]]

DIM_5_TR_SYM_MATRIX = [[9, -2, -2, -2, -2],
                       [-2, 9, -2, -2, -2],
                       [-2, -2, 9, -2, -2],
                       [-2, -2, -2, 9, -2],
                       [-2, -2, -2, -2, 9]]

DIM_7_TR_SYM_MATRIX = [[25, -4, -4, -4, -4, -4, -4],
                       [-4, 25, -4, -4, -4, -4, -4],
                       [-4, -4, 25, -4, -4, -4, -4],
                       [-4, -4, -4, 25, -4, -4, -4],
                       [-4, -4, -4, -4, 25, -4, -4],
                       [-4, -4, -4, -4, -4, 25, -4],
                       [-4, -4, -4, -4, -4, -4, 25]]

DIM_11_TR_SYM_MATRIX = [[21, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2],
                        [-2, 21, -2, -2, -2, -2, -2, -2, -2, -2, -2],
                        [-2, -2, 21, -2, -2, -2, -2, -2, -2, -2, -2],
                        [-2, -2, -2, 21, -2, -2, -2, -2, -2, -2, -2],
                        [-2, -2, -2, -2, 21, -2, -2, -2, -2, -2, -2],
                        [-2, -2, -2, -2, -2, 21, -2, -2, -2, -2, -2],
                        [-2, -2, -2, -2, -2, -2, 21, -2, -2, -2, -2],
                        [-2, -2, -2, -2, -2, -2, -2, 21, -2, -2, -2],
                        [-2, -2, -2, -2, -2, -2, -2, -2, 21, -2, -2],
                        [-2, -2, -2, -2, -2, -2, -2, -2, -2, 21, -2],
                        [-2, -2, -2, -2, -2, -2, -2, -2, -2, -2, 21]]

DIM_13_TR_SYM_MATRIX = [[49, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4],
                        [-4, 49, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4],
                        [-4, -4, 49, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4],
                        [-4, -4, -4, 49, -4, -4, -4, -4, -4, -4, -4, -4, -4],
                        [-4, -4, -4, -4, 49, -4, -4, -4, -4, -4, -4, -4, -4],
                        [-4, -4, -4, -4, -4, 49, -4, -4, -4, -4, -4, -4, -4],
                        [-4, -4, -4, -4, -4, -4, 49, -4, -4, -4, -4, -4, -4],
                        [-4, -4, -4, -4, -4, -4, -4, 49, -4, -4, -4, -4, -4],
                        [-4, -4, -4, -4, -4, -4, -4, -4, 49, -4, -4, -4, -4],
                        [-4, -4, -4, -4, -4, -4, -4, -4, -4, 49, -4, -4, -4],
                        [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 49, -4, -4],
                        [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 49, -4],
                        [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 49]]

# Below is MAGMA code to generate the symmetric matrix of the trace form
"""
p := 3; n := 7; 
R<a0,a1,a2> := PolynomialRing(Integers(),p); 
f := n*(a0^2+a1^2+a2^2)-Integers()!((n-1)/p)*(a0+a1+a2)^2; 
M := SymmetricMatrix(f);
"""
"""
p := 5; n := 11; 
R<a0,a1,a2,a3,a4> := PolynomialRing(Integers(),p); 
f := n*(a0^2+a1^2+a2^2+a3^2+a4^2)-Integers()!((n-1)/p)*(a0+a1+a2+a3+a4)^2; 
M := SymmetricMatrix(f);
"""
"""
p := 7; n := 29; 
R<a0,a1,a2,a3,a4,a5,a6> := PolynomialRing(Integers(),p); 
f := n*(a0^2+a1^2+a2^2+a3^2+a4^2+a5^2+a6^2)-Integers()!((n-1)/p)*(a0+a1+a2+a3+a4+a5+a6)^2; 
M := SymmetricMatrix(f);
"""
"""
p := 11; n := 23; 
R<a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10> := PolynomialRing(Integers(),p); 
f := n*(a0^2+a1^2+a2^2+a3^2+a4^2+a5^2+a6^2+a7^2+a8^2+a9^2+a10^2)-Integers()!((n-1)/p)*(a0+a1+a2+a3+a4+a5+a6+a7+a8+a9+a10)^2; 
M := SymmetricMatrix(f);
"""
"""
p := 13; n := 53; 
R<a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12> := PolynomialRing(Integers(),p); 
f := n*(a0^2+a1^2+a2^2+a3^2+a4^2+a5^2+a6^2+a7^2+a8^2+a9^2+a10^2+a11^2+a12^2)-Integers()!((n-1)/p)*(a0+a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)^2; 
M := SymmetricMatrix(f);
"""

if __name__ == '__main__':
    pass
