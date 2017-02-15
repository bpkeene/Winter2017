set autoscale
unset log
set xtic auto
set ytic auto
set xlabel "Lennard-Jones Time"
set ylabel "<delta r^2>"
set title "Plot of Random Displacement Squared, sim5 msd.dat, dt = 0.005"
plot "sim5_msd.dat" using 2:3 title "Average Squared Displacement" with lines
f(x) = a*x + b
fit f(x) './sim5_msd.dat' u 2:3 via a, b
