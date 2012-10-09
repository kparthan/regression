# Gnuplot script file for plotting data in file "data"

set terminal post eps
set autoscale	# scale axes automatically
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
set title "TRIANGLE"
set xlabel "x"
set ylabel "predictions"
set xr [-0.5:6.78]
set yr [-1.49896:1.49811]
set output "temp/file_XY.eps"
set multiplot
plot "temp/data_XY.txt" using 1:2 title 'original data(y)' \
with points lc rgb "red", \
"temp/data_XY.txt" using 1:3 title 'regression fit(y_estimate)' \
with points lc rgb "blue"
