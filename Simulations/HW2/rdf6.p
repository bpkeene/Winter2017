# Gnuplot script file for plotting the rdf given in Nu_5_step_*0000.rdf
# this file is called rdf5.p
set autoscale
unset log
unset label
set xtic auto
set ytic auto
set title "Radial Distribution Function, nu = 6"
set xlabel "Distance (r/d0)"
set ylabel "g(r)"
set size ratio 1
plot "Nu_6_step_3000000_sim_2.rdf" with lines;

