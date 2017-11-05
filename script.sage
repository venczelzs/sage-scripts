# script.sage SageMath script for solving
# differential equations with Lie symmetry method.

R = PolynomialRing(SR, ['ksi00','ksi10','ksi01','eta00','eta10','eta01'])
P = LaurentPolynomialRing(R, ['x','y','p'])

(ksi00,ksi10,ksi01,eta00,eta10,eta01) = R.gens()
(x,y,p) = P.gens()

# The ODE to be solved is
# y' = f(x,y) = y^2 - y / x (explicit form)
# F(x,y,y') = y' - f(x,y) = y' - y^2 + y / x = 0 (implicit form)
# Let p := y' for convenience.

f = y^2 - y*x^(-1)
F = p - f

ksi = ksi00 + ksi10*x + ksi01*y
eta = eta00 + eta10*x + eta01*y
eta1 = diff(eta,x) + (diff(eta,y) - diff(ksi,x))*p - diff(ksi,y)*p^2

# The infinitesimal generator of the differential operator exp(eps*X).

X = lambda f: ksi*diff(f,x) + eta*diff(f,y) + eta1*diff(f,p)

de = X(F).subs(p=f) # The determining equation.

# From here we want to construct the matrix representing the linear
# system of equations describing X(F) = 0 to solve for ksi and eta.

de_coeffs = de.coefficients()
poly_coeffs = [ksi00,ksi10,ksi01,eta00,eta10,eta01]
m = len(de_coeffs)
n = len(poly_coeffs)

# We construct an m × n matrix whose entry int he i'th row and j'th
# column is the coefficient of the j'th element of poly_coeffs (that
# is the j'th in (ksi00, ..., eta01)) in the i'th element of de_coeffs.

A = matrix(m, n, lambda i, j: QQ(de_coeffs[i].coefficient(poly_coeffs[j])))

print A.rref()
