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
					#J.append(E.isogeny(M[0]**(M[1]-1-r)*(P)).codomain().j_invariant())
					J.append(E.isogeny(M[0]**(M[1]-1-r)*(a*P+b*Q)).codomain().j_invariant())
					#Jb.append(Eb.isogeny(M[0]**(M[1]-1-r)*(IP)).codomain().j_invariant())
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

def test_invariance(l,P,Q,q,MatrixInput):
	'''
	Return the change of basis for which the action of the Frobenius is 
	identical in it

	Input:
	-P,Q a basis of the l**k torsion we study
	-l an integer (which is not useful to this algorithm)
	-q the power of the frobenius
	-MatrixInput the shape of the matrix we want to find thanks to linear
	combinations of P and Q
	
	Output:
	All the change of basis applied to Q->aP + bQ such that we find the 
	same matrix as the input matrix in this new basis.  
	'''
	E=P.curve()
	M=list(Q.order().factor())[0]
	L=[]
	for b in range(Q.order()):
		if (b%M[0]!=0):
			for a in range(P.order()):
				Matrix=etude_action_Frobenius(P,a*P+b*Q,q)
				if(Matrix==MatrixInput):
					print a, " ",b, "\n"; #on imprime les matrices différentes
					L.append([a,b])
					
	return L

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
		TP=k*P #point temporaire d'ordre \ell**r
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

def test_matrices_possibles_sp(L,P,Q,q):
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
				for r in range(2):
					J.append(E.isogeny((M[0]**r)*(P)).codomain().j_invariant())
					J.append(E.isogeny((M[0]**r)*(a*P+b*Q)).codomain().j_invariant())
				#J3=E.isogeny(M[0]**(M[1]-3)*(a*P+b*Q)).codomain().j_invariant()
				if(([Matrix,J] in L)!=True):
					print Matrix, "\n"; #on imprime les matrices différentes
					L.append([Matrix,J])
					L2.append([Matrix,a,b])				
					Test=True
				J=[]
	return [Test,L,L2]

def constructor_isogeny(E,l,q):
	'''
	A function to compute an \ell isogeny starting from E and defined over
	F_{q}
	'''
	L=E(0).division_points(l);
	L=filter(lambda x: x.order()==l, L)
	L=filter(lambda x: ((x[0]^q==x[0]) and (x[1]^q==x[1])), L)
	return E.isogeny(L[randint(0,len(L)-1)])

def retrouve_points_matrice(Ma,q,P,Q):
	'''
	Input:
	-M a matrix that represents the action of the Frobenius
	-q the power of the Frobenius
	-P, Q two basis points of the torsion of the curve

	Ouput:
	-Two points in which the Frobenius acts according to the matrix
	'''	
	E=P.curve()
	M=list(Q.order().factor())[0]
	Mi=matrix([[0,1],[1,0]])*Ma*matrix([[0,1],[1,0]])
	LMa=[]
	LMi=[]
	for b in range(Q.order()):
		if (b%M[0]!=0):
			for a in range(P.order()):
				Matrix=etude_action_Frobenius(P,a*P+b*Q,q)
				if(Matrix==Ma):
					#on teste si les points ne sont pas des multiples d'autres déjà présents
					test=True
					for l in LMa:
						if(l[0].weil_pairing(P,P.order()).multiplicative_order()==1 and l[1].weil_pairing(a*P+b*Q,(a*P+b*Q).order()).multiplicative_order()==1):
							test=False
					#si des multiples ne sont pas deja presents on les ajoute a la liste					
					if(test):
						LMa.append([P,a*P+b*Q])
				if(Matrix==Mi):
					test=True
					for l in LMi:
						if(l[1].weil_pairing(P,P.order()).multiplicative_order()==1 and l[0].weil_pairing(a*P+b*Q,(a*P+b*Q).order()).multiplicative_order()==1):
							test=False
					#si des multiples ne sont pas deja presents on les ajoute a la liste					
					if(test):
						LMi.append([a*P+b*Q,P])
	if(len(LMa)==0 and len(LMi)==0):	
		return False
	else:
		return [LMa,LMi]

def retrouve_points_matrice_improved(Ma,q,P,Q):
	'''
	Input:
	-M a matrix that represents the action of the Frobenius
	-q the power of the Frobenius
	-P, Q two basis points of the torsion of the curve

	Ouput:
	-Two points in which the Frobenius acts according to the matrix
	'''	
	E=P.curve()
	M=list(Q.order().factor())[0]
	Mi=matrix([[0,1],[1,0]])*Ma*matrix([[0,1],[1,0]])
	LMa=[]
	LMi=[]
	for b in range(Q.order()):
		if (b%M[0]!=0):
			for a in range(P.order()):
				Matrix=etude_action_Frobenius(P,a*P+b*Q,q)
				if(Matrix==Ma):
					#on teste si les points ne sont pas des multiples d'autres déjà présents
					test=True
					for l in LMa:
						if(l[0].weil_pairing(P,P.order()).multiplicative_order()==1 and l[1].weil_pairing(a*P+b*Q,(a*P+b*Q).order()).multiplicative_order()==1):
							test=False
					#si des multiples ne sont pas deja presents on les ajoute a la liste					
					if(test):
						LMa.append([P,a*P+b*Q])
				if(Matrix==Mi):
					test=True
					for l in LMi:
						if(l[1].weil_pairing(P,P.order()).multiplicative_order()==1 and l[0].weil_pairing(a*P+b*Q,(a*P+b*Q).order()).multiplicative_order()==1):
							test=False
					#si des multiples ne sont pas deja presents on les ajoute a la liste					
					if(test):
						LMi.append([a*P+b*Q,P])
	P=P+Q
	for b in range(P.order()):
		if (b%M[0]!=0):
			for a in range(Q.order()):
				Matrix=etude_action_Frobenius(Q,a*Q+b*P,q)
				if(Matrix==Ma):
					#on teste si les points ne sont pas des multiples d'autres déjà présents
					test=True
					for l in LMa:
						if(l[0].weil_pairing(Q,Q.order()).multiplicative_order()==1 and l[1].weil_pairing(a*Q+b*P,(a*Q+b*P).order()).multiplicative_order()==1):
							test=False
					#si des multiples ne sont pas deja presents on les ajoute a la liste					
					if(test):
						LMa.append([Q,a*Q+b*P])
				if(Matrix==Mi):
					test=True
					for l in LMi:
						if(l[1].weil_pairing(Q,Q.order()).multiplicative_order()==1 and l[0].weil_pairing(a*Q+b*P,(a*Q+b*P).order()).multiplicative_order()==1):
							test=False
					#si des multiples ne sont pas deja presents on les ajoute a la liste					
					if(test):
						LMi.append([a*Q+b*P,Q])	
	if(len(LMa)==0 and len(LMi)==0):	
		return False
	else:
		return [LMa,LMi]

def test_liste(L,P,Q):
	'''
	Fonction pour tester la sortie de retrouve matrice
	Input:
	-L a list composed of 2 lists of points
	-P,Q two basis points that we are looking for in the list
	Output:
	-La comparaison des points de la liste avec P et Q à l'aide du couplage
	'''
	L0=L[0]
	L1=L[1]
	for l in L0:
		print l[0].weil_pairing(P,P.order()).multiplicative_order(), 'test avec P avec l[0] L0'
		print l[1].weil_pairing(Q,Q.order()).multiplicative_order(), 'test avec Q avec l[1] L0'
		print l[0].weil_pairing(Q,Q.order()).multiplicative_order(), 'test avec Q avec l[0] L0'
		print l[1].weil_pairing(P,P.order()).multiplicative_order(), 'test avec P avec l[1] L0'
		print 'on teste un nouveau couple de points', '\n'
	print '\n', 'on teste la partie symétrique', '\n';
	for l in L1:
		print l[0].weil_pairing(Q,Q.order()).multiplicative_order(), 'test avec Q avec l[0] L1'
		print l[1].weil_pairing(P,P.order()).multiplicative_order(), 'test avec P avec l[1] L1'
		print l[0].weil_pairing(P,P.order()).multiplicative_order(), 'test avec Q avec l[0] L0'
		print l[1].weil_pairing(Q,Q.order()).multiplicative_order(), 'test avec P avec l[1] L0'
		print 'on teste un nouveau couple de points', '\n'

def test_liste_profondeur(L,P,Q,profondeur,b=False):
	'''
	Fonction pour tester la sortie de retrouve matrice
	Input:
	-L a list composed of 2 lists of points
	-P,Q two basis points that we are looking for in the list
	-profondeur an integer that represent the degree of the l-torsion we want to find
	-b a boolean to compute a full test
	'''
	L0=L[0]
	L1=L[1]
	M=list(P.order().factor())[0]
	for l in L0:
		l0=l[0]*M[0]**(M[1]-profondeur)
		l1=l[1]*M[0]**(M[1]-profondeur)
		Pt=P*M[0]**(M[1]-profondeur)
		Qt=Q*M[0]**(M[1]-profondeur)
		print (l0).weil_pairing(Pt,Pt.order()).multiplicative_order(), 'test avec P avec l[0] L0'
		print (l1).weil_pairing(Qt,Qt.order()).multiplicative_order(), 'test avec Q avec l[1] L0'
		if (b):
			print (l0).weil_pairing(Qt,Qt.order()).multiplicative_order(), 'test avec Q avec l[0] L0'
			print (l1).weil_pairing(Pt,Pt.order()).multiplicative_order(), 'test avec P avec l[1] L0'
		print 'on teste un nouveau couple de points', '\n'
	print '\n', 'on teste la partie symétrique', '\n';
	for l in L1:
		l0=l[0]*M[0]**(M[1]-profondeur)
		l1=l[1]*M[0]**(M[1]-profondeur)
		Pt=P*M[0]**(M[1]-profondeur)
		Qt=Q*M[0]**(M[1]-profondeur)
		if (b):
			print (l0).weil_pairing(Qt,Qt.order()).multiplicative_order(), 'test avec Q avec l[0] L1'
			print (l1).weil_pairing(Pt,Pt.order()).multiplicative_order(), 'test avec P avec l[1] L1'
		print (l0).weil_pairing(Pt,Pt.order()).multiplicative_order(), 'test avec P avec l[0] L0'
		print (l1).weil_pairing(Qt,Qt.order()).multiplicative_order(), 'test avec Q avec l[1] L0'
		print 'on teste un nouveau couple de points', '\n'

