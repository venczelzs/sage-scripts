# script.sage SageMath script for solving
# differential equations with Lie symmetry method.

k = 1

# We do this for all k until the matrix A at the bottom has nontrivial kernel,
# that is this should be a while cycle (while (A.kernel().dimension() != 0): k += 1 ...)

var_func = lambda s, k: [(' '.join(('%s%i%i' % (s, i, j) if i + j <= k else '') for j in [0..k] for i in [0..k]))]

igcksi = var_func('ksi', k)
igceta = var_func('eta', k)

inf_gen_coeffs = var(igcksi[0] + igceta[0])
variables = var('x','y','p')

R = PolynomialRing(SR, list(inf_gen_coeffs))
P = LaurentPolynomialRing(R, list(variables))

(ksi00,ksi10,ksi01,eta00,eta10,eta01) = R.gens()
(x,y,p) = P.gens()

# The ODE to be solved is
# y' = f(x,y) = y^2 - y / x (explicit form)
# F(x,y,y') = y' - f(x,y) = y' - y^2 + y / x = 0 (implicit form)
# We denote y' by p for convenience.

f = y^2 - y*x^(-1)
F = p - f

ksi = ksi00 + ksi10*x + ksi01*y
eta = eta00 + eta10*x + eta01*y
eta1 = diff(eta,x) + (diff(eta,y) - diff(ksi,x))*p - diff(ksi,y)*p^2

# The infinitesimal generator of the differential operator exp(eps*X).

X = lambda f: ksi*diff(f,x) + eta*diff(f,y) + eta1*diff(f,p)

determining_equation = X(F).subs(p=f) # The determining equation.

# From here we want to construct the matrix representing the linear
# system of equations describing X(F) = 0 to solve for ksi and eta.

de_coeffs = determining_equation.coefficients()
poly_coeffs = R.gens() # [ksi00,ksi10,ksi01,eta00,eta10,eta01]
m = len(de_coeffs)
n = len(poly_coeffs)

# We construct an m × n matrix whose entry int he i'th row and j'th
# column is the coefficient of the j'th element of poly_coeffs (that
# is the j'# TODO: h in (ksi00, ..., eta01)) in the i'th element of de_coeffs.

A = matrix(m, n, lambda i, j: QQ(de_coeffs[i].coefficient(poly_coeffs[j])))

if A.kernel().dimension() != 0:
    print "nontrivial"
    print A.rref()
else:
    print "trivial"
