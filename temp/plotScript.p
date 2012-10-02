# Gnuplot script file for plotting data in file "data"

set terminal png small
set autoscale	# scale axes automatically
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
set title "SQUARE"
set xlabel "x"
set ylabel "predictions"
set xr [-0.5:6.78]
set yr [-2.17709:2.1935]
set output "temp/file_XY.png"
set multiplot
plot "temp/data_XY.txt" using 1:2 title 'f(x)' \
with points lc rgb "red", \
"temp/data_XY.txt" using 1:3 title 'f(x)+e' \
with points lc rgb "blue"
