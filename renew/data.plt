set size 1.0,1.0
set grid
#set xrange [-2.0:2.0]
#set yrange [-2.0:2.0]
set zrange [-0.1:0.1]
set title "motion"
set hidden3d


splot "position.txt" using 2:3:4 title "motion1" with lines linetype 1 linewidth 2,\
      "position.txt" using 6:7:8 title "motion2" with lines linetype 2 linewidth 2,\
      "position.txt" using 10:11:12 title "motion3" with lines linetype 3 linewidth 2
