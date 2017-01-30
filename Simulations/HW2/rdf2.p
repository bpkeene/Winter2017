# Gnuplot script file for plotting the rdf given in Nu_5_step_100001.rdf
# this file is called rdf5.p
set autoscale
unset log
unset label
set xtic auto
set ytic auto
set title "Radial Distribution Function, nu = 2"
set xlabel "Distance (dimensional)"
set ylabel "g(r)"
set size ratio 1
plot "Nu_2_step_100000.rdf" with lines;

