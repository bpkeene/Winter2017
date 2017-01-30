# Gnuplot script file for plotting the lattice given in *.xyz
# This file is called: latticePicture.p
set autoscale
unset log
unset label
set xtic auto
set ytic auto
set title "Positions of Atoms"
set xlabel "X position"
set ylabel "Y position"
set xr [0.0:1.0]
set yr [0.0:1.0]
set size ratio 1.0
plot "Nu_5_step_100001.xyz" using 1:2 notitle with circles;



