set terminal post eps color enhanced
set autoscale	# scale axes automatically
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
set style data linespoints

set xtics font "Times-Roman, 20"
set ytics font "Times-Roman, 20"
set y2tics font "Times-Roman, 20"

set ytics nomirror
set y2tics nomirror textcolor lt 1;

set xlabel "Number of terms\n"
set ylabel "Message length (in thousands of bits)\n"

set xlabel font "Times-Roman, 25"
set ylabel font "Times-Roman, 25"
set y2label font "Times-Roman, 25"

set key center top font ",20" spacing 1.5 
set key box
set key width 3 

set output "individual_msglens_parabola_sines.eps"

set xrange [1:100]
set yrange [:12]
#set yrange [70:]
#set y2range [50:1000]

set format y '%.2f '
set format y2 '%.2f '

plot "./parabola/sines/1000/results_n1000_s0.2_m1.txt" using 1:($3/1000) title "first part" lc rgb "red" axes x1y2, \
     "" using 1:($4/1000) title "second part" lc rgb "blue" axes x1y1, \
     "" using 1:($5/1000) title "total" lc rgb "dark-green" axes x1y1
