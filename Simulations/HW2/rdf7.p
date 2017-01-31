# Gnuplot script file for plotting the rdf given in Nu_7_step_100001.rdf
# this file is called rdf7.p
set autoscale
unset log
unset label
set xtic auto
set ytic auto
set title "Radial Distribution Function, nu = 7"
set xlabel "Distance"
set ylabel "g(r)"
plot "Nu_7_step_3000000_sim_3.rdf" with lines;

