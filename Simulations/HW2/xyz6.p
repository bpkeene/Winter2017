set autoscale
unset log
unset label
unset key
set xtic auto
set ytic auto
set title "Positions of Atoms - Nu = 6, End of Run"
set xlabel "X position"
set ylabel "Y position"
set xr [0.0:1.0]
set yr [0.0:1.0]
set size ratio 1.0
plot "Nu_6_step_3000000.xyz" w circles ;



