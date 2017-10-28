var('x y p')

F(x,y,p) = p + y

w0 = SR.wild(0)
w1 = SR.wild(1)
w2 = SR.wild(2)

var('ksi00 ksi10 ksi01 eta00 eta10 eta01')

ksi = taylor(function('ksi')(x,y),(x,0),(y,0),1)
ksi = ksi.subs(w0+w1*x+w2*y==ksi00+ksi10*x+ksi01*y)

eta = taylor(function('eta')(x,y),(x,0),(y,0),1)
eta = eta.subs(w0+w1*x+w2*y==eta00+eta10*x+eta01*y)

eta1 = eta.diff(x)+(eta.diff(y)-ksi.diff(x))*p-ksi.diff(y)*p^2

X = lambda f: ksi*f.diff(x)+eta*f.diff(y)+eta1*f.diff(p)

determining_equation = X(F(x,y,p)).subs(p=solve(F(x,y,p),p)[0].rhs())

print determining_equation
print determining_equation.taylor((x,0),(y,0),2)

