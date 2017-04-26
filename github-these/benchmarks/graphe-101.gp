set terminal epslatex color font "Times,16" size 10cm,7cm
#set terminal pdfcairo color enhanced font "Times,16" size 10cm,7cm
set out "graphe-101.tex"

set key inside right top Right
set ylabel "secondes" offset 5
set xlabel "r"
set logscale y 2
set logscale x 2
set xtics 4
plot [4:4096] 'test-script-101-bis.tsv' using 1:($2+$3) title "Algorithmes 1+2" with lines lt 2 dt 2 lw 3 , \
     'test-script-101-bis.tsv' using 1:4  title "Interpolation" with lines lt 1 dt 1 lw 3,
