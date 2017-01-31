# Gnuplot script file for plotting the rdf given in Nu_7_step_100001.rdf
# this file is called rdf7.p
set autoscale
unset log
unset label
set xtic auto
set ytic auto
set title "Power Spectrum, nu = 7"
set xlabel "Index"
set ylabel "Power Spectrum"
plot "Nu_7_step_3000000_sim_3.sf" with lines;

