#!/usr/bin/gnuplot

set terminal pngcairo size 1600,1400 enhanced dashed font 'Verdana,10'
set output "convergence.png"

set multiplot layout 3, 1 title "SaDE pela metado ao 1ZDD"

set style line 11 lc rgb '#808080' lt 1
set border 3 front ls 11
set tics nomirror

set style line 12 lc rgb'#101010' lt 0 lw 1
set grid back ls 12

set style line 1 lw 1 lt 1 lc rgb '#1B9E77' # dark teal
set style line 2 lw 1 lt 1 lc rgb '#D95F02' # dark orange
set style line 3 lw 1 lt 1 lc rgb '#7570B3' # dark lilac
set style line 4 lw 1 lt 1 lc rgb '#E7298A' # dark magenta
set style line 5 lw 1 lt 1 lc rgb '#66A61E' # dark lime green
set style line 6 lw 1 lt 1 lc rgb '#E6AB02' # dark banana
set style line 7 lw 1 lt 1 lc rgb '#A6761D' # dark tan
set style line 8 lw 1 lt 1 lc rgb '#666666' # dark gray

#set title 'Convergence'

set xlabel 'Function Evaluations'
set ylabel 'Score3'

#set key off
set border 3

#set xtics 45

set xrange[0:]
#set yrange[0:1]

plot 'stats.dat' using 2:3 title 'best'   with lines, \
     'stats.dat' using 2:4 title 'mean'   with lines, \


set xlabel 'Function Evaluations'
set ylabel 'Diversity'

plot 'stats.dat' using 2:5 notitle with lines, \
     'stats.dat' using 2:6 notitle with lines, \

set xlabel 'Function Evaluations'
set ylabel 'RMSD'

plot 'stats.dat' using 2:7 notitle with lines, \
