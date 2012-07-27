# Gnuplot script file for plotting data in file "data"

set terminal png small
set autoscale	# scale axes automatically
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
set title "SAWTOOTH"
set xlabel "x"
set ylabel "predictions"
set xr [-1.5:1.5]
set yr [-1.6586:1.65384]
set output "file_XY.png"
set multiplot
plot "data_XY.txt" using 1:2 title 'f(x)' \
with linespoints lc rgb "red", \
"data_XY.txt" using 1:3 title 'f(x)+e' \
with linespoints lc rgb "blue"
