# Gnuplot script file for plotting data in file "data"

set terminal png small
set autoscale	# scale axes automatically
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
set title "SQUARE"
set xlabel "x"
set ylabel "predictions"
set xr [-500.5:500.5]
set yr [-5.10418:5.32109]
set output "file_XY.png"
set multiplot
plot "data_XY.txt" using 1:2 title 'f(x)' \
with linespoints lc rgb "red", \
"data_XY.txt" using 1:3 title 'f(x)+e' \
with linespoints lc rgb "blue"
