set size 1.0,1.0
set grid
set xrange [*:*]
set yrange [*:*]
set title "LPF Effect"



plot   "signal.txt" using 1:2 title "input" with lines linetype 1 linewidth 2,\
       "signal.txt" using 1:3 title "output" with lines linetype 2 linewidth 2,\
	   "signal.txt" using 1:4 title "expectation" with lines linetype 3 linewidth 2

