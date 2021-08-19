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
set key left
set xrange [10:]
set title " CF4 in faujasite, 300K "
plot "0.cf4" using 3:2 title "00%", "10.cf4" using 3:2 title "10%","20.cf4" using 3:2 title "20%", "30.cf4" using 3:2 title "30%","40.cf4" using 3:2 title "40%", "50.cf4" using 3:2 title "50%","60.cf4" using 3:2 title "60%","70.cf4" using 3:2 title "70%","80.cf4" using 3:2 title "80%", "90.cf4" using 3:2 title "90%","100.cf4" using 3:2 title "100%"
pause -1 "Hit Return To Continue"

set title " Methane in faujasite, 300K "
plot "0.me" using 3:2 title "00%", "10.me" using 3:2 title "10%","20.me" using 3:2 title "20%", "30.me" using 3:2 title "30%","40.me" using 3:2 title "40%", "50.me" using 3:2 title "50%","60.me" using 3:2 title "60%","70.me" using 3:2 title "70%","80.me" using 3:2 title "80%", "90.me" using 3:2 title "90%","100.me" using 3:2 title "100%"
pause -1 "Hit Return To Continue"

set terminal postscript landscape color 18

set output "cf4-isos.ps"
set title " CF4 in faujasite, 300K "
plot "0.cf4" using 3:2 title "00%", "10.cf4" using 3:2 title "10%","20.cf4" using 3:2 title "20%", "30.cf4" using 3:2 title "30%","40.cf4" using 3:2 title "40%", "50.cf4" using 3:2 title "50%","60.cf4" using 3:2 title "60%","70.cf4" using 3:2 title "70%","80.cf4" using 3:2 title "80%", "90.cf4" using 3:2 title "90%","100.cf4" using 3:2 title "100%"

set output "meth-isos.ps"
set title " Methane in faujasite, 300K "
plot "0.me" using 3:2 title "00%", "10.me" using 3:2 title "10%","20.me" using 3:2 title "20%", "30.me" using 3:2 title "30%","40.me" using 3:2 title "40%", "50.me" using 3:2 title "50%","60.me" using 3:2 title "60%","70.me" using 3:2 title "70%","80.me" using 3:2 title "80%", "90.me" using 3:2 title "90%","100.me" using 3:2 title "100%"
