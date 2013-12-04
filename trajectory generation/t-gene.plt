set size 1.0,1.0
set grid
set xrange [*:*]
set yrange [*:*]
set title "trajectory generation"



plot   "t-gene.txt" using 1:2 title "input" with lines linetype 1 linewidth 2,\
       "t-gene.txt" using 1:3 title "output" with lines linetype 2 linewidth 2
