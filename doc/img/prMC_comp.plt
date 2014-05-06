set term png font ", 18"
set output "prMC-comp.png"

set xrange [-5:70]

set xtics ("1" 1, "4" 4, "16" 16, "32" 32, "64" 64) font ",18"
set ytics font ",18"

#set xtics scale 2
#set ytics scale 2
#set logscale y


set xlabel "Number of Pthreads" 
set ylabel "Compute Time (s)"

set key right bottom
set border lw 3

plot \
"< sort -nk 3 prMC_N32.dat" u 3:($6>0? $6:1/0) w lp  lw 2 title "32 Nodes" ,\
"< sort -nk 3 prMC_N64.dat" u 3:($6>0? $6:1/0) w lp  lw 2 title "64 Nodes" ,\
"< sort -nk 3 prMC_N128.dat" u 3:($6>0? $6:1/0) w lp lw 2 title "128 Nodes" ,\
"< sort -nk 3 prMC_N256.dat" u 3:($6>0? $6:1/0) w lp lw 2 title "256 Nodes" 

