set terminal x11 enhanced 1
set size 1.0,1.0
set grid
set xrange [-2.0:2.0]
set yrange [-2.0:2.0]
set zrange [-2.0:2.0]
set title "position"
set hidden3d

splot "position.txt" using 2:3:4 title "handPosRef" with lines linetype 1 linewidth 4

set terminal x11 enhanced 2
splot "position.txt" using 5:6:7 title "handPosRes" with lines linetype 2 linewidth 4
