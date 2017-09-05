2017-08-16-presentation_these.pdf: 2017-08-16-presentation_these.tex \
Images/uvsq-logo-cmjn.jpg Images/digiteo.jpg\
Images/graphe-101-149-269-C2.pdf Images/volcan.pdf
	latexmk -pdf 2017-08-16-presentation_these.tex

Images/graphe-101-149-269-C2.pdf: Images/comparaison_C2.tsv \
Images/test-script-18bits-bis.tsv Images/test-script-30bits-bis.tsv \
Images/test-script-62bits-2nd-bis.tsv Images/test-script-252bits-bis.tsv\
Images/graphe-101-149-269-C2.gp
	cd Images && gnuplot graphe-101-149-269-C2.gp


