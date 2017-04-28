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
P=L[1]
M=filter(lambda x: x.weil_pairing(P,2).multiplicative_order()==2,L)
Q=M[1]
for r in range(2):
    L=P.division_points(2)
    P=filter(lambda x: x.order()==2**(r+2),L)[0]
    L=Q.division_points(2)
    Q=filter(lambda x: x.order()==2**(r+2),L)[0]
#Maintenant on a une base de la 8-torsion, plus grande torsion rationnelle sur F_{109**2} et non définie sur F_{109}
#Remarque la 4 torsion est définie sur F_{109}
#Le polynome characterisitique du Frobenius est ******
#Le discriminant vaut ************
#Racines dans CC **************


def etude_action_Frobenius(P,Q,q):
	'''
	Input:
	-P,Q two points that form a basis of the \ell**k torsion such that they
	are not rationnal on the basis field, they are of the same order
	-q the power of the Frobenius

	Output:
	The representation of the action of the Frobenius on this basis as a 
	matrix
	'''
	E=P.curve()
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
	
def test_matrices_possibles(L,P,Q,q,profondeur):
	'''
	Algorithme très gourmand en capacités qui teste toutes les matrices
	possibles pour le Frobenius et compare le résultat à la liste passée en
	paramètre en entrée, si le résultat n'est pas présent à la liste il y est
	rajouté
	Input:
	-L a List of matrices that represent the possible action of the Frobenius
	q on the 2**k-torsion
	-P,Q two points that forms a basis of the 2**k-torsion
	-q the Frobenius we want to test

	Output:
	A boolean value saying that if there is a basis in which the Frobenius
	is different from a matrix from the list and a list concatenated with
	the matrices not present in it 
	'''
	Test=False	
	E=P.curve()
	M=list(Q.order().factor())[0]
	L2=[]
	J=[]
	for b in range(Q.order()):
		if (b%M[0]!=0):
			for a in range(P.order()):
				Matrix=etude_action_Frobenius(P,a*P+b*Q,q)
				for r in range(profondeur):
					J.append(E.isogeny(M[0]**(M[1]-1-r)*(P)).codomain().j_invariant())
					J.append(E.isogeny(M[0]**(M[1]-1-r)*(a*P+b*Q)).codomain().j_invariant())
				#J3=E.isogeny(M[0]**(M[1]-3)*(a*P+b*Q)).codomain().j_invariant()
				if(([Matrix,J] in L)!=True):
					print Matrix, "\n"; #on imprime les matrices différentes
					L.append([Matrix,J])
					L2.append([Matrix,a,b])				
					Test=True
				J=[]
	return [Test,L,L2]
#Pour travailler avec la 16-torsion
F1.<a>=FiniteField(149^4)
E1=EllipticCurve(j=F1(56))
if(E1.cardinality()!=492860160):
	E1=E1.quadratic_twist()
L=filter(lambda x: x.order()==2, E1(0).division_points(2))
P1=L[1]
M=filter(lambda x: x.weil_pairing(P1,2).multiplicative_order()==2,L)
Q1=M[1]
for r in range(3):
    L=P1.division_points(2)
    P1=filter(lambda x: x.order()==2**(r+2),L)[0]
    L=Q1.division_points(2)
    Q1=filter(lambda x: x.order()==2**(r+2),L)[0]			
#R.<c> = Zq(4, prec = 20); 
#Contruction de l'extension 2-adic non ramifiée à l'aide de Sage
def constructor_isogeny(E,l,q):
	'''
	A function to compute an \ell isogeny starting from E and defined over
	F_{q}
	'''
	L=E(0).division_points(l);
	L=filter(lambda x: x.order()==l, L)
	L=filter(lambda x: ((x[0]^q==x[0]) and (x[1]^q==x[1])), L)
	return E.isogeny(L[0])

def test_ultime(l,P,Q,q,profondeur):
	'''
	A function that test if the mapping of the isogeny preserves the action
	of the Frobenius

	Input:
	-P,Q two basis points of a curve
	-q the Frobenius to test
	-l the degree of the isogeny we want to test
	-profondeur the depth we want to test in the volcano
	'''
	Test=False
	E=P.curve()
	M=list(Q.order().factor())[0]
	I=constructor_isogeny(E,l,q)
	L2=[]
	IP=I(P)
	IQ=I(Q)
	L=[]
	Eb=IP.curve()
	J=[]
	Jb=[]
	for b in range(Q.order()):
		if (b%M[0]!=0):
			for a in range(P.order()):
				Matrix=etude_action_Frobenius(P,a*P+b*Q,q)
				Matrixb=etude_action_Frobenius(IP,a*IP+b*IQ,q)
				for r in range(profondeur):
					J.append(E.isogeny(M[0]**(M[1]-1-r)*(P)).codomain().j_invariant())
					J.append(E.isogeny(M[0]**(M[1]-1-r)*(a*P+b*Q)).codomain().j_invariant())
					Jb.append(Eb.isogeny(M[0]**(M[1]-1-r)*(IP)).codomain().j_invariant())
					Jb.append(Eb.isogeny(M[0]**(M[1]-1-r)*(a*IP+b*IQ)).codomain().j_invariant())
				#J1=E.isogeny(M[0]**(M[1]-1)*(a*P+b*Q)).codomain().j_invariant()
				#J2=E.isogeny(M[0]**(M[1]-2)*(a*P+b*Q)).codomain().j_invariant()
				#J3=E.isogeny(M[0]**(M[1]-1)*(a*P+b*Q)).codomain().j_invariant()
				#J4=E.isogeny(M[0]**(M[1]-2)*(a*P+b*Q)).codomain().j_invariant()
				#J1b=Eb.isogeny(M[0]**(M[1]-1)*(a*IP+b*IQ)).codomain().j_invariant()
				#J2b=Eb.isogeny(M[0]**(M[1]-2)*(a*IP+b*IQ)).codomain().j_invariant()
				#J3=E.isogeny(M[0]**(M[1]-3)*(a*P+b*Q)).codomain().j_invariant()
				if(([Matrix,J,Matrixb,Jb] in L)!=True):
					print Matrix, "\n"; #on imprime les matrices différentes
					L.append([Matrix,J,Matrixb,Jb])
					L2.append([Matrix,a,b])				
					Test=True
				J=[]
				Jb=[]
	return [Test,L,L2]

