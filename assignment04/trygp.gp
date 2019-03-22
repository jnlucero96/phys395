#!/usr/bin/gnuplot
#
# Generate an animated spiral
#
# AUTHOR: Hagen Wierstorf

reset

# png
set terminal pngcairo size 350,262 enhanced font 'Verdana,10'

# color definitions
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2 # --- blue

unset key
# set border 0
# unset tics
# set view 342,0
set xrange [0:100]
set xlabel "time"
set yrange [-4:4]
set ylabel "pend. 2 pos."
system('mkdir -p png')
# spiral upwards
n=0
do for [ii=1:10000:100] {
    n=n+1
    set output sprintf('png/spiral%03.0f.png',n)
    plot 'DATA' u 1:3 every ::1::ii w l ls 1, \
         'DATA' u 1:3 every ::ii::ii w p ls 1
}