#code pour une courbe située sur le cratère d'un volcan de Atkin, le volcan a la forme suivante avec les j-invariant:
#	le volcan est régulier sur F_{149^2}	
#
#		     56	
#      /	     |		   \	
#    7		    35            114	
#/	\	/	\	/	\
#81	112	123	110	101	141
#les courbes sonf définies sur F_149	
F.<a>=FiniteField(149^2)
E=EllipticCurve(j=F(56))
if(E.cardinality()!=22464):
	E=E.quadratic_twist()
L=filter(lambda x: x.order()==2, E(0).division_points(2))
P=L[randint(0,len(L)-1)]
M=filter(lambda x: x.weil_pairing(P,2).multiplicative_order()==2,L)
Q=M[randint(0,len(M)-1)]
for r in range(2):
    L=P.division_points(2)
    P=filter(lambda x: x.order()==2**(r+2),L)[randint(0,len(L)-1)]
    L=Q.division_points(2)
    Q=filter(lambda x: x.order()==2**(r+2),L)[randint(0,len(L)-1)]
#Maintenant on a une base de la 8-torsion, plus grande torsion rationnelle sur F_{109**2} et non définie sur F_{109}
#Remarque la 4 torsion est définie sur F_{109}
#Le polynome characterisitique du Frobenius est ******
#Le discriminant vaut ************
#Racines dans CC **************


#Pour travailler avec la 16-torsion
F1.<a>=FiniteField(149^4)
E1=EllipticCurve(j=F1(56))
if(E1.cardinality()!=492860160):
	E1=E1.quadratic_twist()
L=filter(lambda x: x.order()==2, E1(0).division_points(2))
P1=L[randint(0,len(L)-1)]
M=filter(lambda x: x.weil_pairing(P1,2).multiplicative_order()==2,L)
Q1=M[randint(0,len(M)-1)]
for r in range(3):
    L=P1.division_points(2)
    P1=filter(lambda x: x.order()==2**(r+2),L)[randint(0,len(L)-1)]
    L=Q1.division_points(2)
    Q1=filter(lambda x: x.order()==2**(r+2),L)[randint(0,len(L)-1)]			
#R.<c> = Zq(4, prec = 20); 
#Contruction de l'extension 2-adic non ramifiée à l'aide de Sage


#Pour travailler avec la 32 torsion
F2.<a>=FiniteField(149^8)
E2=EllipticCurve(j=F2(56))
if(E2.cardinality()!=242935033147223040):
	E2=E2.quadratic_twist()
L=filter(lambda x: x.order()==2, E2(0).division_points(2))
P2=L[randint(0,len(L)-1)]
M=filter(lambda x: x.weil_pairing(P2,2).multiplicative_order()==2,L)
Q2=M[randint(0,len(M)-1)]
for r in range(4):
    	L=P2.division_points(2)
    	P2=filter(lambda x: x.order()==2**(r+2),L)[randint(0,len(L)-1)]
    	L=Q2.division_points(2)
  	Q2=filter(lambda x: x.order()==2**(r+2),L)[randint(0,len(L)-1)]	

#Pour travailler sur une courbe plus bas dans le volcan
Lb=E1(0).division_points(2);
Lb=filter(lambda x: x.order()==2, Lb)
Eb=E1.isogeny(Lb[-1]).codomain()
Lb=filter(lambda x: x.order()==2 ,Eb(0).division_points(2))
Pb=Lb[randint(0,len(Lb)-1)]
Lb=filter(lambda x: x.weil_pairing(Pb,2).multiplicative_order()==2, Lb)
Qb=Lb[randint(0,len(Lb)-1)]
for r in range(2):
	Lb=Pb.division_points(2)
	Pb=filter(lambda x: x.order()==2**(r+2),Lb)[randint(0,len(Lb)-1)]
	Lb=Qb.division_points(2)
	Qb=filter(lambda x: x.order()==2**(r+2),Lb)[randint(0,len(Lb)-1)]

#Test pour l iosgenie avec le Frobenius
Ma=etude_action_Frobenius(P,Q,149)
print Ma
#On calcule la matrice associée à la base dont on cherche un candidat pour l'image
phi=constructor_isogeny(P.curve(),3,149)
print phi, phi.codomain().j_invariant()
Pc=phi(P)
Qc=phi(Q)
Pcache=Pc
Qcache=Qc
Pc,Qc=Qc,Pc
Qc=3*Qc+4*Pc
print Qc.weil_pairing(Pc,Pc.order()).multiplicative_order()==Pc.order()
Listest=retrouve_points_matrice(Ma,149,Pc,Qc)
Listest2=retrouve_points_matrice(Ma,149,Qc,Pc)
Listest3=retrouve_points_matrice(Ma,149,Qc,Pc+Qc)

