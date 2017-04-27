#code pour une courbe située sur le cratère d'un volcan de Atkin, le volcan a la forme suivante avec les j-invariant:
#		
#
#		     45	
#      /	     |		   \	
#    1		   0            54	
#les courbes sonf définies sur F_109	
F.<a>=FiniteField(109^2)
E=EllipticCurve(j=F(45))
E.cardinality().factor()

