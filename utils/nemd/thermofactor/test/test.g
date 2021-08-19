#set linestyle 1 lt 1 lw 3 pt 1 ps 20
#set linestyle 2 lt 2 lw 3 pt 2 ps 20
#   set linestyle 10 lt 6 lw 6
#   plot "isotherm.Methane" with l ls 10
#   set pointsize 1.5
# f(x)=k*x/(c-x)
# k=741
# c=18
# fit f(x) 'isotherm.Methane' using 2:1 via k, c
# plot "isotherm.Methane" using 1:2 title "simulaton",  c*x/(k+x) title "Langmuir Fit"
#plot "isotherm.n_Heptane" using 1:2 with l ls 1, "igcmc.dat" using 1:2 with p ps 2




set data style linespoints

set xlabel "Total Pressure, fug. kPa "
set ylabel "molec/ uc "
set logscale x

set title " CF4 in faujasite, 300K "
plot "test.dat" using 1:2 title "sim","test.dat" using 1:3 title "fit"
pause -1 "Hit Return To Continue"
