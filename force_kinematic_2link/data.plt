set size 1.0,1.0
set grid
set xrange [-2.0:2.0]
set yrange [-2.0:2.0]
set zrange [-2.0:2.0]
set title "motion"
set hidden3d

set parametric
set angle degree
set urange [0:360]
set vrange [0:360]
set isosample 36,36
set ticslevel 0
set size 0.7,1.0
a=1.0

splot 'force_res.dat' using 5:6:7 title "motion" with lines linetype 1 linewidth 2,\
      a*cos(u)*cos(v),a*sin(u)*cos(v),a*sin(v) notitle linestyle 5 linewidth 1
