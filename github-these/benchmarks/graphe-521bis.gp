#set terminal epslatex color font "Times,16" size 20cm,7cm
#set terminal pdfcairo color enhanced font "Times,16" size 10cm,7cm
#set out "graphe-521.tex"
set terminal pdf color enhanced font "Times,16"
set out "graphe-521bis.pdf"

set key inside right bottom Right
set ylabel "secondes" offset 5
set xlabel "r"
set logscale y 2
set logscale x 2
set xtics 4
plot [4:4096] 'test-script-521-bis.tsv' using 1:4 with line dt 1 lw 3 title "p = 521", \
     'test-script-1033-bis.tsv' using 1:4 with line dt 4 lw 3 title "p = 1033", \
#     'test-script-62bits-bis.tsv' using 1:4 with line dt 6 lw 3 title "p ≈ 2^{62}", \
#     'test-script-30bits-bis.tsv' using 1:4 with line dt 2 lw 3 title "$p \simeq 2^{30}$", \
#     'test-script-62bits-2nd-bis.tsv' using 1:4 with line dt 4 lw 3 title "$p \simeq 2^{62}$", \
#     'test-script-252bits-bis.tsv' using 1:4 with line dt 5 lw 3 title "$p \simeq 2^{252}$"

#     'test-script-149-bis.tsv' using 1:4 with line lt 2 lw 3 title "p=149", \
#     'test-script-269-bis.tsv' using 1:4 with line lt 3 lw 3 title "p=269", \
#     'test-script-521-bis.tsv' using 1:4 with line lt 4 lw 3 title "p=521", \
#     'test-script-1033-bis.tsv' using 1:4 with line lt 5 lw 3 title "p=1033", \
