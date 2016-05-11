load('../isogeny-computing-testless.sage')
load('../tate-module-extension-testless.sage')
load('../extension_corps.sage')

n=15
r=10
Fichier = open('../../benchmarks/test-script-252bits-bis.tsv','a')

p=7237005577332262213973186563042994240829374041602535252466099000494570602917
F=FiniteField(p)
K=Tower_two(F,1)
print K._top
E=EllipticCurve(j=F(72))
if E.cardinality()!=7237005577332262213973186563042994240736408655068204765837505660044652071200:
	E=E.quadratic_twist()
print E
for end in range(2,443):
	deg=end**2
	a,b,c,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	a,b2,c2,d,e,f,g=tate_module(E,((16.0)/3)*deg,K,2,conservation=True)
	if calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)[0]==True:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init_bis(b,c,b,c,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init_bis(b,c,b,c,2,d,deg,e,f,g)',number=n,repeat=r,seconds=True,preparse=True)
		A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=n,repeat=r,seconds=True,preparse=True)		
	else:
		M,Lc,P2,Q2,R2,o2,o1,Tower,deg=calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)
		C=timeit('calcul_isogenie_init(b,c,b,c,2,d,deg,e,f,g)',number=n,repeat=r,seconds=True,preparse=True)		
		if calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)[0]==True:
			A=timeit('calcul_isogenie_step_bis(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=n,repeat=r,seconds=True,preparse=True)
		else:
			A=timeit('calcul_isogenie_step(M,Lc,P2,Q2,R2,o2,o1,Tower,deg)',number=n,repeat=r,seconds=True,preparse=True)
	B=timeit('tate_module(E,((16.0)/3)*deg,K,2,conservation=True)',number=n,repeat=r,seconds=True,preparse=True)
	Fichier.write( str(deg) +'\t'+ str(B) +'\t'+ str(C) +'\t'+ str(A) + '\t'+ str(d) +'\n')
	Fichier.close()

