set autoscale
unset log
set xtic auto
set ytic auto
set xlabel "Step (number)"
set ylabel "Energy"
set title "Plot of Energy, sim6\_prop.dat, dt = 0.005"
plot "sim6_prop.dat" using 1:2 title "Potential Energy" with lines, \
"sim6_prop.dat" using 1:3 title "Kinetic Energy" with lines, \
"sim6_prop.dat" using 1:4 title "Total Energy" with lines
