set size 1.0,1.0
set grid
set xrange [-1.0:1.0]
set yrange [-1.0:1.0]
set zrange [-1.0:1.0]
set title "COM position"
set hidden3d


splot "COMpos.txt" using 2:3:4 title "posRef" with lines linetype 1 linewidth 3,\
      "COMpos.txt" using 5:6:7 title "posRes" with lines linetype 2 linewidth 3
