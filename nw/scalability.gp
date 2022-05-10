set term svg size 1920,1080
set size ratio 0.8
set output 'scalability.svg'
set title "Strong scaling test per version of the Needleman-Wunsch algorithm"
set title font ', 35'

set key left top Right

# set key at 2,0.9,0

ideal(x)=x
one(x)=1

set key box lw 2 vertical
set key samplen 5 height 0.6 autotitle
set key spacing 1 font ',25'

set xlabel 'Number of processes' font 'Helvetica Bold, 30' offset 0,-1
set xtics font ', 30' offset 0,-0.5

set ylabel 'Speedup' font 'Helvetica Bold, 30' offset -0.5
set ytics font ', 30'

# set autoscale xfix
# set autoscale yfix
# set ytics 2

set logscale x 2
set logscale y 2

# set yrange [0.5:*]

# set xrange [1:*]

# set offsets 0, 0.1, 0.05, 0 # (left, right, top, bottom)
# set offsets 0.1, 0.1

set style line 1 lt 1 lw 4 ps 3 pt 1


# plot 'plot_height.dat' using ($3/$2):xticlabel(1) smooth bezier with lines title 'Average efficiency' ls 1, \
#                    '' using ($3/$2):xticlabel(1) lw 50 ps 2.3 title 'Real Efficiency', \
#                    1 title "Ideal Efficiency"

# '' with labels center offset 3.4,.5 notitle



plot \
      ideal(x) title "Ideal speedup" ls 1 lc 4, \
      one(x) title "bin's speedup" ls 1 lc 2, \
      one(x) title "char's speedup" ls 1 lc 2, \
      'scalability.dat' using 1:($3/$2) with lp ls 1 lc 3 title "diag's speedup"