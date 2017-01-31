# Gnuplot script file for plotting the rdf given in Nu_7_step_100001.rdf
# this file is called rdf7.p
set autoscale
unset log
unset label
set xtic auto
set ytic auto
set title "Free Energy Surface, nu = 7"
set xlabel "Distance"
set ylabel "F/kT"
plot "Nu_7_step_3000000_sim_3.fes" with lines;

