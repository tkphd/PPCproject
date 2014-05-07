set term png font ", 18"
set output "pr-speed.png"

#set xtics ("1" 1, "4" 4, "16" 16, "32" 32, "64" 64) font ",18"
set ytics font ",18"

#set xtics scale 2
#set ytics scale 2
#set logscale x
#set logscale y

set xlabel "MPI ranks"
set ylabel "Speedup (arb. units)"

set key right bottom
set border lw 3

plot \
"< awk '$3==1' pr.dat | sort -nk 1" u ($1):(1/($8)) w lp  lw 2 title " 1 thread" ,\
"< awk '$3==4' pr.dat | sort -nk 1" u ($1):(1/($8)) w lp  lw 2 title " 4 pthreads" ,\
"< awk '$3==16' pr.dat | sort -nk 1" u ($1):(1/($8)) w lp  lw 2 title "16 pthreads" ,\
"< awk '$3==32' pr.dat | sort -nk 1" u ($1):(1/($8)) w lp  lw 2 title "32 pthreads" ,\
"< awk '$3==64' pr.dat | sort -nk 1" u ($1):(1/($8)) w lp  lw 2 title "64 pthreads"

