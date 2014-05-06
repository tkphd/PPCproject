set term png font ", 18"
set output "prMC-initBW.png"

set xrange [-5:70]
#set yrange [0:16]

set xtics ("1" 1, "4" 4, "16" 16, "32" 32, "64" 64) font ",18"
set ytics font ",18"

#set xtics scale 2
#set ytics scale 2
#set ytics 2 
#set logscale y


set xlabel "Number of Pthreads" 
set ylabel "Init Bandwidth (MB/s)"

set border lw 3

set key right top

f(x) = x/1048576

plot \
"< sort -nk 3 prMC_N32.dat" u 3:(f($5)>0? f($5): 1/0) w lp  lw 2 title "32 Nodes" ,\
"< sort -nk 3 prMC_N64.dat" u 3:(f($5)>0? f($5): 1/0) w lp  lw 2 title "64 Nodes" ,\
"< sort -nk 3 prMC_N128.dat" u 3:(f($5)>0? f($5): 1/0) w lp lw 2 title "128 Nodes" ,\
"< sort -nk 3 prMC_N256.dat" u 3:(f($5)>0? f($5): 1/0) w lp lw 2 title "256 Nodes" 

