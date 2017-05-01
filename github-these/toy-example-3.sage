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
#Dans F3
F3.<a>=FiniteField(149^16)
E1b=EllipticCurve(j=F3(56))
if E1b.cardinality()!=59017430136820283374698267657154560:
	E1b=E1b.quadratic_twist()
Lb2=E1b(0).division_points(2)
L2b=filter(lambda x: x.order()==2, L2b)
E2b=E1b.isogeny(Lb[-1]).codomain()
L2b=filter(lambda x: x.order()==2 ,E2b(0).division_points(2))
P2b=L2b[randint(0,len(L2b)-1)]
L2b=filter(lambda x: x.weil_pairing(P2b,2).multiplicative_order()==2, L2b)
Q2b=L2b[randint(0,len(L2b)-1)]
for r in range(3):
	L2b=P2b.division_points(2)
	P2b=filter(lambda x: x.order()==2**(r+2),L2b)[randint(0,len(L2b)-1)]
	L2b=Q2b.division_points(2)
	Q2b=filter(lambda x: x.order()==2**(r+2),L2b)[randint(0,len(L2b)-1)]

#Dans F2
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

#Test pour l isogenie avec le Frobenius sur la 8 torsion
Ma=etude_action_Frobenius(P,Q,149)
#On calcule la matrice associée à la base dont on cherche un candidat pour l'image
phi=constructor_isogeny(P.curve(),3,149)
print phi, phi.codomain().j_invariant()
Pc=phi(P)
Qc=phi(Q)
Pcache=Pc
Qcache=Qc
Pc,Qc=Qc,Pc #on permute l ordre de la base pour brouiller les pistes
Qc=3*Qc+2*Pc
print Qc.weil_pairing(Pc,Pc.order()).multiplicative_order()==Pc.order()
Listest=retrouve_points_matrice(Ma,149,Pc,Qc)
Listest2=retrouve_points_matrice(Ma,149,Qc,Pc)
Listest3=retrouve_points_matrice(Ma,149,Qc,Pc+Qc)
Listest4=retrouve_points_matrice(Ma,149,Pc+Qc,Qc)

#Listesti=retrouve_points_matrice_improved(Ma,149,Pc,Qc)
#Listesti2=retrouve_points_matrice_improved(Ma,149,Qc,Pc)
#Listesti3=retrouve_points_matrice_improved(Ma,149,Qc,Pc+Qc)
#Listesti4=retrouve_points_matrice_improved(Ma,149,Pc+Qc,Qc)


#Test pour l isogenie avec le Frobenius sur la 16 torsion
Ma1=etude_action_Frobenius(P1,Q1,149)
#On calcule la matrice associée à la base dont on cherche un candidat pour l'image
phi1=constructor_isogeny(P1.curve(),3,149)
print phi1, phi1.codomain().j_invariant()
Pc1=phi1(P1)
Qc1=phi1(Q1)
Pcache1=Pc1
Qcache1=Qc1
Pc1,Qc1=Qc1,Pc1 #on permute l ordre de la base pour brouiller les pistes
Qc1=3*Qc1+2*Pc1
print Qc1.weil_pairing(Pc1,Pc1.order()).multiplicative_order()==Pc1.order()
Listest1=retrouve_points_matrice(Ma,149,Pc1,Qc1)
Listest12=retrouve_points_matrice(Ma,149,Qc1,Pc1)
Listest13=retrouve_points_matrice(Ma,149,Qc1,Pc1+Qc1)
Listest14=retrouve_points_matrice(Ma,149,Pc1+Qc1,Qc1)

#Listesti1=retrouve_points_matrice_improved(Ma,149,Pc1,Qc1)
#Listesti12=retrouve_points_matrice_improved(Ma,149,Qc1,Pc1)
#Listesti13=retrouve_points_matrice_improved(Ma,149,Qc1,Pc1+Qc1)
#Listesti14=retrouve_points_matrice_improved(Ma,149,Pc1+Qc1,Qc1)

#Test avec la connaissance d'une corresepondance de sous groupe
#L=(2*Pcache1).division_points(2)
Pc1=phi1(P1)
Qc1=phi1(Q1)
L=(2*(Pc1+2*Qc1)).division_points(2)
L1=filter(lambda x: x.weil_pairing(Pcache1,x.order()).multiplicative_order()==1, L)
L2=filter(lambda x: x.weil_pairing(Pcache1,x.order()).multiplicative_order()!=1, L)

Pc1=phi1(P1)
Qc1=phi1(Q1)

Qc1=3*Qc1+2*Pc1
if(len(L1)>1):
	Listestame1=retrouve_points_matrice(Ma,149,L1[0],Qc1)
	Listestame12=retrouve_points_matrice(Ma,149,L1[1],Qc1)
	Listestame13=retrouve_points_matrice(Ma,149,Qc1,L1[0])
	Listestame14=retrouve_points_matrice(Ma,149,Qc1,L1[1])

if(len(L2)>1):
	Listestameb1=retrouve_points_matrice(Ma,149,L2[0],Qc1)
	Listestameb12=retrouve_points_matrice(Ma,149,L2[1],Qc1)
	Listestameb13=retrouve_points_matrice(Ma,149,Qc1,L2[0])
	Listestameb14=retrouve_points_matrice(Ma,149,Qc1,L2[1])

