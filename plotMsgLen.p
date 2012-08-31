set term post eps
set autoscale	# scale axes automatically
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
set title "N = 1000, Sigma = 0.75"
set xlabel "# of terms"
set ylabel "Message Length"
set output "./results/results_n1000_s0.75.txt.eps"
plot "./results/results_n1000_s0.75.txt" using 1:3 notitle with linespoints lc rgb "blue"
