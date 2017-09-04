2016-these.pdf: 2016-these.tex Biblio.bib \
benchmarks/graphe-101bis.pdf benchmarks/graphe-101-2^18bis.pdf\
benchmarks/graphe-521bis.pdf modele_these_cyril.pdf
	makeindex 2016-these.nlo -s nomencl.ist -o 2016-these.nls
	latexmk -pdf 2016-these.tex

modele_these_cyril.pdf: modele_these_cyril.tex 
	latexmk -pdf modele_these_cyril.tex

benchmarks/graphe-101bis.pdf: benchmarks/test-script-101-bis.tsv benchmarks/graphe-101bis.gp
	cd benchmarks && gnuplot graphe-101bis.gp

benchmarks/graphe-101-2^18bis.pdf: benchmarks/test-script-101-bis.tsv \
benchmarks/test-script-18bits-bis.tsv benchmarks/test-script-30bits-bis.tsv \
benchmarks/test-script-62bits-2nd-bis.tsv benchmarks/test-script-252bits-bis.tsv\
benchmarks/graphe-101-2^18bis.gp
	cd benchmarks && gnuplot graphe-101-2^18bis.gp

benchmarks/graphe-521bis.pdf: benchmarks/test-script-521-bis.tsv \
benchmarks/test-script-1033-bis.tsv benchmarks/graphe-521bis.gp
	cd benchmarks && gnuplot graphe-521bis.gp
