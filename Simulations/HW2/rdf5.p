# Gnuplot script file for plotting the rdf given in Nu_5_step_*0000.rdf
# this file is called rdf5.p
set autoscale
unset log
unset label
set xtic auto
set ytic auto
set title "Radial Distribution Function, nu = 5"
set xlabel "Distance (r/d0)"
set ylabel "g(r)"
set size ratio 1
plot "Nu_5_step_20000.rdf" with lines;

