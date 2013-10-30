set size 1.0,1.0
set grid
set xrange [-2.0:2.0]
set yrange [-2.0:2.0]
set title "compos"



plot "COMpos.txt" using 2:3 title "compos" with lines linetype 1 linewidth 5
