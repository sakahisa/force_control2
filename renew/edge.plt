set size 1.0,1.0
set grid
set xrange [-3.0:3.0]
set yrange [-3.0:3.0]
set title "edgepos"
set size square


plot   "edgepos.txt" using 2:3 title "COMpos" with lines linetype 1 linewidth 5,\
       "edgepos.txt" using 7:8 title "edgepos" with lines linetype 2 linewidth 5
