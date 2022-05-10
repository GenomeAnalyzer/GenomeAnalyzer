set term svg size 1920,1080
set size ratio 0.8
set output 'cycles.svg'
set title "Number of execution cycles per sequence size and version of the Needleman-Wunsch algorithm"
set title font 'Helvetica, 35'

set key left top Right

# set key at 2,0.9,0

set key box lw 2 vertical
set key samplen 5 height 0.6 autotitle
set key spacing 1 font ',25'

set xlabel "Sequence size" font 'Helvetica Bold, 35' offset 0,-3
set xtics font ', 30' offset 0,-0.5

set ylabel "Number of execution cycles" font 'Helvetica Bold, 35' offset -10,-1.5
set ytics font ', 30'

set autoscale xfix
set autoscale yfix
# set yrange [0:1.05]
# set offsets 0.1, 0.1, 0.05, 0 # (left, right, top, bottom)
set offsets 0.1, 0.1

set style line 1 lt 1 lw 4 ps 3 pt 1

# plot 'plot_height.dat' using ($3/$2):xticlabel(1) smooth bezier with lines title 'Average efficiency' ls 1, \
#                    '' using ($3/$2):xticlabel(1) lw 50 ps 2.3 title 'Real Efficiency', \
#                    1 title "Ideal Efficiency"

# '' with labels center offset 3.4,.5 notitle


plot 'bin.out' using ($1*12*16):3:4 with lp ls 1 lc 1 title "bin version", \
    'char.out' using ($1*12*16):3:4 with lp ls 1 lc 2 title "char version", \
    'diag.out' using ($1*12*16):3:4 with lp ls 1 lc 3 title "diag version"