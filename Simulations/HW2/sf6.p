# Gnuplot script file for plotting the rdf given in Nu_7_step_100001.rdf
# this file is called rdf7.p
set autoscale
unset label
set xtic auto
set ytic auto
set title "Structure Factor, nu = 6"
set xlabel "Inverse Distance (1/ k*d0)"
set ylabel "Structure Factor"
plot "Nu_6_step_3000000_sim_2.sf" with lines;

