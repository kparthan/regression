# Gnuplot script file for plotting data in file "Results"

set terminal png small
set autoscale	# scale axes automatically
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
#set title "SAWTOOTH"
set xlabel "M"
set ylabel "Error/Message Length"
set xr [1:100]
#set yr [-4.80207:5.12866]
set output "results.png"
set multiplot
plot "results_1000.txt" using 2:4 title 'error' \
with linespoints lc rgb "red"
plot "results_1000.txt" using 2:5 title 'msgLen' \
with linespoints lc rgb "blue"
