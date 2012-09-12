# Gnuplot script file for plotting data in file "data"

set terminal png small
set autoscale	# scale axes automatically
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
set title "SAWTOOTH"
set xlabel "x"
set ylabel "f(x)"
set xr [-1.5:1.5]
<<<<<<< HEAD:test/plotScript.p
set yr [-0.5:1.5]
set output "file_XfX.png"
=======
set yr [-1.5:1.5]
set output "file_XY.png"
>>>>>>> 57d61fedf3014331915abae9f24ac44ddf0b747b:plotScript.p
set multiplot
plot "fun.txt" using 1:2 title 'f(x)' \
with linespoints lc rgb "red"

