set size 1.0,1.0
set grid
set xrange [*:*]
set yrange [*:*]
set zrange [*:*]
set title "COM position"
set hidden3d


splot "position.txt" using 2:3:4 title "posRef" with lines linetype 1 linewidth 3,\
      "position.txt" using 5:6:7 title "posRes" with lines linetype 2 linewidth 3,\
      "position.txt" using 8:9:10 title "COM_start1" with lines linetype 3 linewidth 3,\
      "position.txt" using 12:13:14 title "COM_start2" with lines linetype 4 linewidth 3,\
      "position.txt" using 16:17:18 title "COM_start3" with lines linetype 5 linewidth 3
