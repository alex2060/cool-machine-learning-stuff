import sympy
from sympy import *
from sympy import I
from sympy import sympify
# put c(1,2,3) in for s(1,2,3) and r(1,2,3) in for r(1,2,3)
# i as I 
h1=1
h2=2
h3=4
c1=1
c2=3
c3=4


x = Symbol('x')
b = Symbol('b')

s1 = Symbol('s1')
s2 = Symbol('s2')
s3 = Symbol('s3')
r1 = Symbol('r1')
r2 = Symbol('r2')
r3 = Symbol('r3')

A=h1/c1
B=h2/c2
C=h3/c3
Vec=[A,B,C]
M=Matrix([[(exp(-2)-1)/2,(exp(-3)-1)/3,(exp(-4)-1)/4],[(exp(-3)-1)/3,(exp(-4)-1)/4,(exp(-5)-1)/5],[(exp(-4)-1)/4,(exp(-5)-1)/5,(exp(-6)-1)/6]])
Minv=M.inv()
#print(Minv)
soltion=Minv.dot(Vec)
#print(sol[1])


m1=r1*exp(I*(x+b))
m2=r2*exp(I*(x+b)*2)
m3=r3*exp(I*(x+b)*3)
m=m1+m2+m3
#print("m=m1+m2+m3=")
#pprint(m1+m2+m3)
#print()
#print()
g=h1*exp(I*b)+h2*exp(I*b*2)+h3*exp(I*b*3)
#print("g=")
#pprint(g)
#print()
#print()
f1=I*s1*exp(I*x)
f2=I*s2*exp(I*x*2)
f3=I*s3*exp(I*x*3)
f=f1+f2+f3
#print("f=f1+f2+f3")
#pprint(f1*(m1))
# intigrating F(x)g(x+1) form 1 to i
end=I
k12=integrate(f1*(m1), (x,0,end))
k14=integrate(f2*(m1), (x,0,end))
k16=integrate(f3*(m1), (x,0,end))

k22=integrate(f1*(m2), (x,0,end))
k24=integrate(f2*(m2), (x,0,end))
k26=integrate(f3*(m2), (x,0,end))

k32=integrate(f1*(m3), (x,0,end))
k34=integrate(f2*(m3), (x,0,end))
k36=integrate(f3*(m3), (x,0,end))

e1=k36+k34+k32
e2=k26+k24+k22
e3=k16+k14+k12


#e3=factor(e3)
print("given g(b)")
pprint(g)
print("given m(x+b)")
pprint(m)
print("matrix to be inversed")
pprint(M)
print("aproximation of marix")
pprint(M.evalf())


# intagral of F(x)g(x+1) form 1 to i with b(1,2,3) and r(1,2,3) left in (unevlauted)
print("intagral of F(x)g(x+1) form 1 to i with s(1,2,3) and r(1,2,3) left in (unevlauted) =") 
print(e3)
print("+")
print(e2)
print("+")
print(e1)
print()
print()


# intagral of F(x)g(x+1) form 1 to i with b(1,2,3) and r(1,2,3) left in aproximated=
print("intagral of F(x)g(x+b) form 1 to i with s(1,2,3) and r(1,2,3) left in aproximated=")
print(factor(e3.evalf()))
print(factor(e2.evalf()))
print(factor(e1.evalf()))
print()
print()


print("intagral of F(x)g(x+b) form 1 to i with aproximated=")
print(factor(e3.subs(r1, c1).subs(s1, soltion[0]).subs(s2, soltion[1]).subs(s3, soltion[2]).evalf()))
print(factor(e2.subs(r2, c2).subs(s1, soltion[0]).subs(s2, soltion[1]).subs(s3, soltion[2]).evalf()))
print(factor(e1.subs(r3, c3).subs(s1, soltion[0]).subs(s2, soltion[1]).subs(s3, soltion[2]).evalf()))
print()
print()

print("F(x) aproximated")
faprox=f.subs(s1, soltion[0]).subs(s2, soltion[1]).subs(s3, soltion[2]).evalf()
pprint(faprox)

maprox=m.evalf().subs(r1, c1).subs(r2, c2).subs(r3, c3)
print("m")
print(maprox)
print("m*f")
mf=sympify(maprox*f.evalf())
print(mf)
print("aproximation of intagral of F(x)g(x+b) form 1 to i-g(x)")
aprox=integrate(faprox*m1.subs(r1,c1), (x,0,end))+integrate(faprox*m2.subs(r2,c2), (x,0,end))+integrate(faprox*m3.subs(r3,c3), (x,0,end))-g
print(aprox.evalf())
input()

