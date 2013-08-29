set size 1.0,1.0
set grid
set xrange [-2.0:2.0]
set yrange [-2.0:2.0]

set title "X position"

set terminal x11 enhanced 1
plot "position.txt" using 1:2 title "posRef" with lines linetype 1 linewidth 4
replot "position.txt" using 1:5 title "posRes" with lines linetype 2 linewidth 4

set terminal x11 enhanced 2
plot "position.txt" using 1:3 title "posRef" with lines linetype 1 linewidth 4
replot "position.txt" using 1:6 title "posRes" with lines linetype 2 linewidth 4

set terminal x11 enhanced 3
plot "position.txt" using 1:4 title "posRef" with lines linetype 1 linewidth 4
replot "position.txt" using 1:7 title "posRes" with lines linetype 2 linewidth 4
