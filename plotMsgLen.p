set term post eps
set autoscale	# scale axes automatically
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
set title "N = 100, Sigma = 0"
set xlabel "# of terms"
set ylabel "Message Length"
set output "./Results/results_n100_s0.txt.eps"
plot "./Results/results_n100_s0.txt" using 1:3 notitle with linespoints lc rgb "blue"
