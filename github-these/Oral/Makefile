2017-08-16-presentation_these.pdf: 2017-08-16-presentation_these.tex \
Images/uvsq-logo-cmjn.jpg Images/digiteo.jpg\
Images/graphe-101-149-269-C2.pdf Images/volcan.pdf\
Images/graphe-101-149-269.pdf
	latexmk -pdf 2017-08-16-presentation_these.tex

Images/graphe-101-149-269-C2.pdf: Images/comparaison_C2-100bits.tsv \
Images/test-script-18bits-bis.tsv Images/test-script-30bits-bis.tsv \
Images/test-script-62bits-2nd-bis.tsv Images/test-script-252bits-bis.tsv\
Images/graphe-101-149-269-C2.gp Images/comparaison_C2-1000bits.tsv
	cd Images && gnuplot graphe-101-149-269-C2.gp

Images/graphe-101-149-269.pdf: Images/graphe-101-149-269.gp \
Images/test-script-18bits-bis.tsv Images/test-script-30bits-bis.tsv \
Images/test-script-62bits-2nd-bis.tsv Images/test-script-252bits-bis.tsv\

	cd Images && gnuplot graphe-101-149-269.gp
