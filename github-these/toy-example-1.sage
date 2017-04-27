#code pour une courbe située sur le cratère d'un volcan de Atkin, le volcan a la forme suivante avec les j-invariant:
#	le volcan est régulier sur F_{109^2}	
#
#		     6	
#      /	     |		   \	
#    65		    22             94	
#/	\	/	\	/	\
#69	106	80	61	85	26
#les courbes sonf définies sur F_109	
F.<a>=FiniteField(109^2)
E=EllipticCurve(j=F(6))
E.cardinality().factor()
L=E(0).division_points(2); L
P=L[1]
M=filter(lambda x: x.weil_pairing(P,2).multiplicative_order()==2,L)
Q=M[1]
for r in range(2):
    L=P.division_points(2)
    P=filter(lambda x: x.order()==2**(r+2),L)[0]
    L=Q.division_points(2)
    Q=filter(lambda x: x.order()==2**(r+2),L)[0]
#Maintenant on a une base de la 8-torsion, plus grande torsion rationnelle sur F_{109**2} et non définie sur F_{109}
#Remarque meme la 4 torsion n'est pas définie sur F_{109}
#Le polynome characterisitique du Frobenius est X**2+214-109**2
#Le discriminant vaut -1728

def etude_action_Frobenius(P,Q,q):
	'''
	Input:
	-P,Q two points that form a basis of the \ell**k torsion such that they
	are not rationnal on the basis field, they are of the same order

	Output:
	The representation of the action of the Frobenius on this basis as a 
	matrix
	'''
	PF=E(P[0]**q,P[1]**q)
	QF=E(Q[0]**q,Q[1]**q)
	L=list(P.order().factor())[0]
	a1=0; 
	b1=0;
	a2=0;
	b2=0;
	#coefficients of the vector representing the action of the Frobenius on P
	k=P.order()
	for r in range(L[1]):
		k=int(k/L[0])
		TP=k*P
		TPF=k*PF
		TQ=k*Q
		TQF=k*QF
		for a0 in range(L[0]):
			for b0 in range(L[0]):
				if(TPF==(a1+a0*L[0]**r)*TP+(b1+b0*L[0]**r)*TQ):
					a1=a1+a0*L[0]**r
					b1=b1+b0*L[0]**r
				if(TQF==(a2+a0*L[0]**r)*TP+(b2+b0*L[0]**r)*TQ):
					a2=a2+a0*L[0]**r
					b2=b2+b0*L[0]**r
	return matrix([[a1,a2],[b1,b2]])
	
	


