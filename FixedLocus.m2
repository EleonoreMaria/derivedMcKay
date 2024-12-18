--G12

k = QQ[t]/ideal(t^8-t^4+1) --t is a 24th root of unity
e = t^3
i = t^6
w = t^8
r = (e+e^7)/2 --1/sqrt(2)
R = k[x,y]
sigma = map(R,R,{i*x,-i*y})
tau = map(R,R,{-y,x})
mu = map(R,R,{r*(e^7*x+e^5*y),r*(e^7*x+e*y)})
alpha = map(R,R,{e*x,e^3*y})
p1 = x^2-y^2
p2 = x^2 + y^2
p3 = x*y
q1 = x^3+(2*w+1)*x*y^2
q2 = y^3+(2*w+1)*x^2*y
s1 = x^3+(2*w^2+1)*x*y^2
s2 = y^3+(2*w^2+1)*x^2*y
gamma1 = x^5-5*x*y^4
gamma2 = y^5-5*x^4*y
phi = p2^2+4*w*p3^2
psi = p2^2+4*w^2*p3^2
T=p1*p2*p3
W=phi*psi

N = ideal(p1*p2*p3,phi*psi,phi^3)
V4rho1' = ideal(phi)
V8rho1' = ideal(psi^2)
V4rho1'' = ideal(psi)
V8rho1'' = ideal(phi^2)
V5rho2' = ideal(x*phi,y*phi)
V7rho2' = ideal(s1*psi,s2*psi)
V5rho2'' = ideal(x*psi,y*psi)
V7rho2'' = ideal(q1*phi,q2*phi)
V5rho2 = ideal(gamma1,gamma2)
V7rho2 = ideal(s1*phi, s2*phi)
V2rho3phi = ideal(x^2*phi,x*y*phi,y^2*phi)
V2rho3psi = ideal(x^2*psi,x*y*psi,y^2*psi)

--verifying H-invariant submodules V_5(\rho_2)\oplus V_7(\rho_2) for E(\rho_2)
tau(gamma1)==-gamma2
tau(s2*phi)==s1*phi
tau(gamma2)==gamma1
tau(s1*phi)==-s2*phi
sigma(gamma1)==i*gamma1
sigma(s2*phi)==i*s2*phi
sigma(gamma2)==-i*gamma2
sigma(s1*phi)==-i*s1*phi
mu(gamma1)==((-i+1)/(2))*gamma1-((i+1)/(2))*gamma2
mu(s2*phi)==((-i+1)/(2))*s2*phi+((i+1)/(2))*s1*phi
mu(gamma2)==((-i+1)/(2))*gamma1+((i+1)/(2))*gamma2
mu(s1*phi)==-((-i+1)/(2))*s2*phi+((i+1)/(2))*s1*phi

--verifying E(\rho_2) fixed pointwise under alpha
I=ideal(gamma1,gamma2,s1*phi, s2*phi)
alpha(I)==I
X= ideal(gamma1+s2*phi,gamma2-s1*phi)
alpha(X)==X
Y= ideal(gamma1-s2*phi,gamma2+s1*phi)
alpha(Y)==Y
Z= ideal(gamma1,gamma2)
alpha(Z)==Z

--finding the intersection of E(\rho_2) and E(\rho_3) and verifying it is fixed under alpha
J=ideal(y*gamma1+x*gamma2,x*gamma1+y*gamma2,x*gamma1-y*gamma2)
tau(J)==J
sigma(J)==J
mu(J)==J
K=ideal(y*gamma1+x*gamma2,x*gamma1+y*gamma2,x*gamma1-y*gamma2,s1*phi, s2*phi)
alpha(K)==K



--E(\rho_3) fixed under alpha, finding H-invariant submodules
L=ideal(x^2*phi,x*y*phi,y^2*phi,x^2*psi,x*y*psi,y^2*psi)
alpha(L)==L
tau(phi)==phi
tau(psi)==psi
tau(x^2+y^2)==x^2+y^2
tau(x*y)==-x*y
tau(x^2-y^2)==-(x^2-y^2)
sigma(phi)==phi
sigma(psi)==psi
sigma(x^2+y^2)==-(x^2+y^2)
sigma(x*y)==-x*y
sigma(x^2-y^2)==-(x^2-y^2)
mu(x^2+y^2)==-i*(x^2-y^2)
mu(x^2-y^2)==-2*x*y
mu(x*y)==-(i/2)*(x^2+y^2)
mu(phi)==w^2*phi
mu(psi)==w*psi

M=ideal((x^2-y^2)*(phi+w*psi),x*y*(w*phi+w*psi),(x^2+y^2)*(w^2*phi+w*psi))
M==J
O=ideal(p1*(-w^2*phi+psi),p3*(-phi+psi),p2*(-w*phi+psi))
alpha(O)==O

--Checking (1,e) is fixed in the non-exceptional locus
alpha(T)==T
alpha(W)==W
alpha(phi)==-psi
sub(T,{x=>1, y=>e})
sub(W,{x=>1, y=>e})
sub(phi,{x=>1, y=>e})
sub(psi,{x=>1, y=>e})
P=ideal(T-2*e,phi^3-(2*t^6-4*t^2)^3,W+12)
Q=ideal(T-2*e,psi^3-(-2*t^6+4*t^2)^3,W+12)
P==Q

--checking intersection of non-exceptional and exceptional locus
restart
s=QQ[a]
S=frac s
k = S[t]/ideal(t^8-t^4+1)
e = t^3 --eighth root of unity
i = t^6 --fourth root of unity
w = t^8 --cube root of unity
r = (e+e^7)/2 --1/sqrt(2)
R = k[x,y]
p1 = x^2-y^2
p2 = x^2 + y^2
p3 = x*y
q1 = x^3+(2*w+1)*x*y^2
q2 = y^3+(2*w+1)*x^2*y
s1 = x^3+(2*w^2+1)*x*y^2
s2 = y^3+(2*w^2+1)*x^2*y
gamma1 = x^5-5*x*y^4
gamma2 = y^5-5*x^4*y
phi = p2^2+4*w*p3^2
psi = p2^2+4*w^2*p3^2
f1 = p1*p2*p3
f2 = phi*psi
f3 = x^12-33*x^8*y^4-33*x^4*y^8+y^12
It=ideal(f1-2*e*a^6,f2+12*a^8,f3)

A = p1*(-w^2*phi+psi)
B = p3*(-phi+psi)
C = p2*(-w*phi+psi)

A%It --observe all terms in remainder are divisible by a
B%It --observe all terms in remainder are divisible by a
C%It --observe all terms in remainder are divisible by a

--G13

k = QQ[t]/ideal(t^8-t^4+1) --adjoins a primitive 24th root of unity to our field
e = t^3 -- primitive 8th root of unity
i = t^6 -- primitive 4th root of unity
w = t^8 -- primitive cube root of unity
r = (e+e^7)/2 -- 1/sqrt(2)
R = k[x,y]
sigma = map(R,R,{i*x,-i*y})
tau = map(R,R,{-y,x})
mu = map(R,R,{r*(e^7*x+e^5*y),r*(e^7*x+e*y)})
kappa = map(R,R,{e*x,e^7*y})
alpha = map(R,R,{e*x,e^3*y})

-- Ito-Nakamura p.227
p1 = x^2-y^2
p2 = x^2 + y^2
p3 = x*y
phi = p2^2+4*w*p3^2
psi = p2^2+4*w^2*p3^2
T = p1*p2*p3
W = phi*psi

-- Ito-Nakamura p.237
X = x^12-33*x^8*y^4-33*x^4*y^8+y^12
F = X*T
N = ideal(W,T^2,F)

-- Ito-Nakamura p.235
-- E(rho1')
V6rho1' = ideal(T)
V12rho1' = ideal(x^12-33*x^8*y^4-33*x^4*y^8+y^12)
alpha(V6rho1'+V12rho1')==V6rho1'+V12rho1'
alpha(V6rho1')==V6rho1'
V7rho2' = ideal(x*T,y*T)
alpha(V12rho1'+V7rho2')==V12rho1'+V7rho2'
Xrho1' = ideal(T+X)
sigma(Xrho1')==Xrho1'
tau(Xrho1')==Xrho1'
mu(Xrho1')==Xrho1'
kappa(Xrho1')==Xrho1'
alpha(Xrho1')==Xrho1'

-- E(rho2)
V7rho2 = ideal(7*x^4*y^3+y^7,-x^7-7*x^3*y^4)
V11rho2 = ideal(x^10*y-6*x^6*y^5+5*x^2*y^9,-x*y^10+6*x^5*y^6-5*x^9*y^2)
alpha(V7rho2+V11rho2)==V7rho2+V11rho2
alpha(V7rho2)==V7rho2
V8rho3 = ideal(-2*x*y^7-14*x^5*y^3,x^8-y^8,2*x^7*y+14*x^3*y^5)
alpha(V11rho2+V8rho3)==V11rho2+V8rho3

f1 = 7*x^4*y^3+y^7
f2 = -x^7-7*x^3*y^4
g1 = x^10*y-6*x^6*y^5+5*x^2*y^9
g2 = -x*y^10+6*x^5*y^6-5*x^9*y^2
Xrho2 = ideal(f1+g1,f2+g2)
sigma(Xrho2)==Xrho2
tau(Xrho2)==Xrho2
mu(Xrho2)==Xrho2
kappa(Xrho2)==Xrho2
alpha(Xrho2)==Xrho2

-- E(rho2')
V11rho2' = ideal(-11*x^8*y^3-22*x^4*y^7+y^11,11*x^3*y^8+22*x^7*y^4-x^11)
alpha(V7rho2'+V11rho2')==V7rho2'+V11rho2'
alpha(V12rho1'+V7rho2')==V12rho1'+V7rho2'
V8rho3' = ideal(x^2*T,x*y*T,y^2*T)
alpha(V11rho2'+V8rho3')==V11rho2'+V8rho3'

v1 = -11*x^8*y^3-22*x^4*y^7+y^11
v2 = 11*x^3*y^8+22*x^7*y^4-x^11
Xrho2' = ideal(x*T+v1,y*T+v2)
sigma(Xrho2)==Xrho2
tau(Xrho2)==Xrho2
mu(Xrho2)==Xrho2
kappa(Xrho2)==Xrho2
alpha(Xrho2)==Xrho2

-- E(rho2'')
V8rho2'' = ideal(psi^2,-phi^2)
V10rho2'' = ideal(x^5*y*psi-x*y^5*phi,-x^5*y*phi+x*y^5*psi)
alpha(V8rho2''+V10rho2'')==V8rho2''+V10rho2''
alpha(V8rho2'')==V8rho2''
S1V8rho2'' = ideal(x*psi^2,y*psi^2,-x*phi^2,-y*phi^2)
alpha(V10rho2''+S1V8rho2'')==V10rho2''+S1V8rho2''

Xrho2'' = ideal(psi^2+T*phi,phi^2+T*psi)
sigma(Xrho2'')==Xrho2''
tau(Xrho2'')==Xrho2''
mu(Xrho2'')==Xrho2''
kappa(Xrho2'')==Xrho2''
alpha(Xrho2'')==Xrho2''

-- E(rho3)
V10rho3 = ideal(4*x^10+60*x^6*y^4,5*x^9*y+54*x^5*y^5+5*x*y^9,60*x^4*y^6+4*y^10)
alpha(V8rho3+V10rho3)==V8rho3+V10rho3
alpha(V11rho2+V8rho3)==V11rho2+V8rho3

h1 = -2*x*y^7-14*x^5*y^3
h2 = x^8-y^8
h3 = 2*x^7*y+14*x^3*y^5
S1V8rho3=ideal(x*h1,x*h2,x*h3,y*h1,y*h2,y*h3)
alpha(V10rho3+S1V8rho3)==V10rho3+S1V8rho3

j1 = 4*x^10+60*x^6*y^4
j2 = 5*x^9*y+54*x^5*y^5+5*x*y^9
j3 = 60*x^4*y^6+4*y^10
Xrho3 = ideal(h1+j1,h2+j2,h3+j3)
sigma(Xrho3)==Xrho3
tau(Xrho3)==Xrho3
mu(Xrho3)==Xrho3
kappa(Xrho3)==Xrho3
alpha(Xrho3)==Xrho3

-- E(rho3')
V8rho3' = ideal(x^2*T,x*y*T,y^2*T)
V10rho3' = ideal(-3*x^8*y^2-14*x^4*y^6+y^10,8*x^7*y^3+8*x^3*y^7,x^10-14*x^6*y^4-3*x^2*y^8)
alpha(V8rho3'+V10rho3')==V8rho3'+V10rho3'
alpha(V11rho2'+V8rho3')==V11rho2'+V8rho3'
S1V8rho3'=ideal(x^3*T,x^2*y*T,x*y^2*T,y^3*T)
alpha(V10rho3'+S1V8rho3')==V10rho3'+S1V8rho3'

u1 = -3*x^8*y^2-14*x^4*y^6+y^10
u2 = 8*x^7*y^3+8*x^3*y^7
u3 = x^10-14*x^6*y^4-3*x^2*y^8
Xrho3' = ideal(x^2*T+u1,x*y*T+u2,y^2*T+u3)
sigma(Xrho3')==Xrho3'
tau(Xrho3')==Xrho3'
mu(Xrho3')==Xrho3'
kappa(Xrho3')==Xrho3'
alpha(Xrho3')==Xrho3'


--G22

k = QQ[s]/ideal(s^8-s^6+s^4-s^2+1) --20th cyclotomic extension of the rationals
e=s^4 --fifth root of unity
t=s^4 --fifth root of unity
i=s^5 --fourth root of unity
r = (2*(e+e^4)+1)/5 -- 1/sqrt(5), note sqrt(5) is 5*r
R = k[x,y]
sigma = map(R,R,{-e^3*x,-e^2*y})
tau = map(R,R,{r*(-(e-e^4)*x+(e^2-e^3)*y),r*((e^2-e^3)*x+(e-e^4)*y)})
alpha = map(R,R,{i*r*(-e+e^4)*x+i*r*(e^2-e^3)*y,i*r*(e^2-e^3)*x+i*r*(e-e^4)*y})

--verification group is G_22
tauinv=inverse(tau)
alpha(tauinv(x))==i*x
alpha(tauinv(y))==i*y

--Ito-Nakamura p.239
s1 = x^10+66*x^5*y^5-11*y^10 --sigma1
s2 = -11*x^10-66*x^5*y^5+y^10 --sigma2
t1 = x^10-39*x^5*y^5-26*y^10 --tau1
t2 = -26*x^10+39*x^5*y^5+y^10 --tau2

--Ito-Nakamura p.241
--E(rho2)
V11rho2 = ideal(x*s1,-y*s2)
V19rho2 = ideal(-57*x^15*y^4+247*x^10*y^9+171*x^5*y^14+y^19,
    	    	-x^19+171*x^14*y^5-247*x^9*y^10-57*x^4*y^15)
c1=x*s1  --from V11rho2
c2=-y*s2  --from V11rho2
d1=-57*x^15*y^4+247*x^10*y^9+171*x^5*y^14+y^19 --from V19rho2
d2=-x^19+171*x^14*y^5-247*x^9*y^10-57*x^4*y^15 --from V19rho2
alpha(V11rho2+V19rho2)==V11rho2+V19rho2
Xrho2 =ideal(c1+d1,c2-d2)
sigma(Xrho2)==Xrho2
tau(Xrho2)==Xrho2
alpha(Xrho2)==Xrho2

--E(rho2')
V13rho2' = ideal(y^3*t2,
                 -x^3*t1)
V17rho2' = ideal(x^17+119*x^12*y^5+187*x^7*y^10+17*x^2*y^15,
                 -17*x^15*y^2+187*x^10*y^7-119*x^5*y^12+y^17)
u1=y^3*t2 -- from V13rho2' 
u2=-x^3*t1 -- from V13rho2'
v1=x^17+119*x^12*y^5+187*x^7*y^10+17*x^2*y^15 -- from V17rho2'
v2=-17*x^15*y^2+187*x^10*y^7-119*x^5*y^12+y^17 -- from V17rho2'                 
alpha(V13rho2'+V17rho2')==V13rho2'+V17rho2'
Xrho2' =ideal(u1+v1,u2+v2)
sigma(Xrho2')==Xrho2'
tau(Xrho2')==Xrho2'
alpha(Xrho2')==Xrho2'

--E(rho3)
V12rho3 = ideal(x^2*s1,
    	     	-5*x^11*y-5*x*y^11,
                y^2*s2)
V18rho3 = ideal(-12*x^15*y^3+117*x^10*y^8+126*x^5*y^13+y^18,
                45*x^14*y^4-130*x^9*y^9-45*x^4*y^14,
		x^18-126*x^13*y^5+117*x^8*y^10+12*x^3*y^15)
f1=x^2*s1  -- from V12rho3
f2=-5*x^11*y-5*x*y^11 -- from V12rho3
f3=y^2*s2 -- from V12rho3
g1=-12*x^15*y^3+117*x^10*y^8+126*x^5*y^13+y^18 --from V18rho3
g2=45*x^14*y^4-130*x^9*y^9-45*x^4*y^14  -- from V18rho3
g3=x^18-126*x^13*y^5+117*x^8*y^10+12*x^3*y^15 -- from V18rho3 
alpha(V12rho3+V18rho3)==V12rho3+V18rho3
Xrho3=ideal(f1+g1,f2+g2,f3+g3)
sigma(Xrho3)==Xrho3
tau(Xrho3)==Xrho3
alpha(Xrho3)==Xrho3
	    
--E(rho3'')	    
V14rho3'' = ideal(x^14-14*x^9*y^5+49*x^4*y^10,
    	    	  7*x^12*y^2-48*x^7*y^7-7*x^2*y^12,
		  49*x^10*y^4+14*x^5*y^9+y^14)
V16rho3'' = ideal(3*x^15*y-143*x^10*y^6-39*x^5*y^11+y^16,
    	    	  -25*x^13*y^3-25*x^3*y^13,
		  x^16+39*x^11*y^5-143*x^6*y^10-3*x*y^15)

m1=x^14-14*x^9*y^5+49*x^4*y^10 -- from v14rho3''
m2=7*x^12*y^2-48*x^7*y^7-7*x^2*y^12 -- from v14rho3''
m3=49*x^10*y^4+14*x^5*y^9+y^14 -- from v14rho3''
n1=3*x^15*y-143*x^10*y^6-39*x^5*y^11+y^16 --from V16rho3''
n2=-25*x^13*y^3-25*x^3*y^13 --from V16rho3''
n3=x^16+39*x^11*y^5-143*x^6*y^10-3*x*y^15 --from V16rho3''

alpha(V14rho3''+V16rho3'')==V14rho3''+V16rho3''
alpha(V14rho3'')==V14rho3''
Xrho3''=ideal(m1+n1,m2+n2,m3+n3)
sigma(Xrho3'')==Xrho3''
tau(Xrho3'')==Xrho3''
alpha(Xrho3'')==Xrho3''

--E(rho4)
V13rho4 = ideal(x^3*s1,
                -3*x^12*y+22*x^7*y^6-7*x^2*y^11,
		-7*x^11*y^2-22*x^6*y^7-3*x*y^12,
		y^3*s2)
V17rho4 = ideal(-2*x^15*y^2+52*x^10*y^7+91*x^5*y^12+y^17,
                10*x^14*y^3-65*x^9*y^8-35*x^4*y^13,
		-35*x^13*y^4+65*x^8*y^9+10*x^3*y^14,
		-x^17+91*x^12*y^5-52*x^7*y^10-2*x^2*y^15)
h1=x^3*s1 --from v13rho4
h2=-3*x^12*y+22*x^7*y^6-7*x^2*y^11 --from v13rho4
h3=-7*x^11*y^2-22*x^6*y^7-3*x*y^12--from v13rho4
h4=y^3*s2--from v13rho4
j1=-2*x^15*y^2+52*x^10*y^7+91*x^5*y^12+y^17 --from v17rho4
j2=10*x^14*y^3-65*x^9*y^8-35*x^4*y^13 --from v17rho4
j3=-35*x^13*y^4+65*x^8*y^9+10*x^3*y^14 --from v17rho4
j4=-x^17+91*x^12*y^5-52*x^7*y^10-2*x^2*y^15 --from v17rho4

alpha(V13rho4+V17rho4)==V13rho4+V17rho4
Xrho4=ideal(h1+j1,h2+j2,h3+j3,h4+j4)
sigma(Xrho4)==Xrho4
tau(Xrho4)==Xrho4
alpha(Xrho4)==Xrho4

--E(rho4')
V14rho4' = ideal(x*y^3*t2,
                 -x^4*t1,
		 y^4*t2,
		 -x^3*y*t1)
V16rho4' = ideal(-2*x^15*y+77*x^10*y^6-84*x^5*y^11+y^16,
    	    	 35*x^12*y^4+110*x^7*y^9+15*x^2*y^14,
		 15*x^14*y^2-110*x^9*y^7+35*x^4*y^12,
		 -x^16-84*x^11*y^5-77*x^6*y^10-2*x*y^15)
a1=x*y^3*t2 --from V14rho4'
a2=-x^4*t1  --from V14rho4'
a3=y^4*t2  --from V14rho4'
a4=-x^3*y*t1  --from V14rho4'
b1=-2*x^15*y+77*x^10*y^6-84*x^5*y^11+y^16  --from V16rho4'
b2=35*x^12*y^4+110*x^7*y^9+15*x^2*y^14  --from V16rho4'
b3=15*x^14*y^2-110*x^9*y^7+35*x^4*y^12  --from V16rho4'
b4=-x^16-84*x^11*y^5-77*x^6*y^10-2*x*y^15  --from V16rho4'	     
	     
alpha(V14rho4'+V16rho4')==V14rho4'+V16rho4'

Xrho4'=ideal(a1+b2,a2+b1,a3+b4,a4+b3)
sigma(Xrho4')==Xrho4'
tau(Xrho4')==Xrho4'
alpha(Xrho4')==Xrho4'

--E(rho5)
V14rho5 = ideal(x^4*s1,
                -2*x^13*y+33*x^8*y^6-8*x^3*y^11,
		-5*x^12*y^2-5*x^2*y^12,
		-8*x^11*y^3-33*x^6*y^8-2*x*y^13,
		-y^4*s2)
V16rho5 = ideal(64*x^15*y+728*x^10*y^6+y^16,
                66*x^14*y^2+676*x^9*y^7-91*x^4*y^12,
		56*x^13*y^3+741*x^8*y^8-56*x^3*y^13,
		91*x^12*y^4+676*x^7*y^9-66*x^2*y^14,
		x^16+728*x^6*y^10-64*x*y^15)
k1=x^4*s1 --from V14rho5
k2=-2*x^13*y+33*x^8*y^6-8*x^3*y^11 --from V14rho5
k3=-5*x^12*y^2-5*x^2*y^12 --from V14rho5
k4=-8*x^11*y^3-33*x^6*y^8-2*x*y^13 --from V14rho5
k5=-y^4*s2 --from V14rho5
l1=64*x^15*y+728*x^10*y^6+y^16 -- from v16rho5
l2=66*x^14*y^2+676*x^9*y^7-91*x^4*y^12 -- from v16rho5
l3=56*x^13*y^3+741*x^8*y^8-56*x^3*y^13 -- from v16rho5
l4=91*x^12*y^4+676*x^7*y^9-66*x^2*y^14 -- from v16rho5
l5=x^16+728*x^6*y^10-64*x*y^15 -- from v16rho5

alpha(V14rho5+V16rho5)==V14rho5+V16rho5
Xrho5=ideal(k1+l1,k2+l2,k3+l3,k4+l4,k5-l5)
sigma(Xrho5)==Xrho5
tau(Xrho5)==Xrho5
alpha(Xrho5)==Xrho5

--checking intersection of non-exceptional and exceptional locus
restart
w=QQ[a]
S=frac w
k = S[s]/ideal(s^8-s^6+s^4-s^2+1) --20th cyclotomic extension of the rationals
e=s^4 --fifth root of unity
t=s^4 --fifth root of unity
i=s^5 --fourth root of unity
r = (2*(e+e^4)+1)/5 -- 1/sqrt(5), note sqrt(5) is 5*r
R = k[x,y]

s1 = x^10+66*x^5*y^5-11*y^10 --sigma1
s2 = -11*x^10-66*x^5*y^5+y^10 --sigma2
t1 = x^10-39*x^5*y^5-26*y^10 --tau1
t2 = -26*x^10+39*x^5*y^5+y^10 --tau2
f1=x*y*(x^10+11*x^5*y^5-y^10)
f2=x^20-228*x^15*y^5+494*x^10*y^10+228*x^5*y^15+y^20
f3=x^30+522*x^25*y^5-10005*x^20*y^10-10005*x^10*y^20-522*x^5*y^25+y^30
It=ideal(f1-s*(2+11*i)*a^12,f2+(492+456*i)*a^20,f3)
A=x^14-14*x^9*y^5+49*x^4*y^10
B=7*x^12*y^2-48*x^7*y^7-7*x^2*y^12
C=49*x^10*y^4+14*x^5*y^9+y^14

A%It --observe all terms in remainder are divisible by a
B%It --observe all terms in remainder are divisible by a
C%It --observe all terms in remainder are divisible by a
