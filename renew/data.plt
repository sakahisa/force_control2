set size 1.0,1.0
set grid
set xrange [*:*]
set yrange [*:*]
set zrange [*:*]
set title "force"
set hidden3d


splot "force.txt" using 2:3:4 title "forceRef" with lines linetype 1 linewidth 5,\
      "force.txt" using 8:9:10 title "forceRes" with lines linetype 2 linewidth 2

