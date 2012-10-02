# Gnuplot script file for plotting data in file "temp/test_msglen"

set terminal post eps
set key left top
set autoscale	# scale axes automatically
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
set title "COMPARISON OF MESSAGE LENGTHS\nN = 100, Sigma = 0"
set xlabel "# of terms"
set ylabel "message length"
#set xr [-1.5:1.5]
#set yr [-0.829709:1.77094]
set output "temp/comparison_msglen.eps"
set multiplot
plot "temp/test_msglen" using 1:6 title 'PART 1' \
with linespoints lc rgb "red", \
"temp/test_msglen" using 1:7 title 'PART 2' \
with linespoints lc rgb "blue", \
"temp/test_msglen" using 1:8 title 'TOTAL' \
with linespoints lc rgb "green"
