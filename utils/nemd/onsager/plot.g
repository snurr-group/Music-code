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

#set logscale x
#set yrange [1:14]

set title " FC_11"
set xlabel "time, ps"
set ylabel "<Vi(0).Vj(t)>, (ang/ps)^2"
plot "FC_11" using 1:2 
pause -1 "Hit Return To Continue"

set title " FC_12"
set xlabel "time, ps"
set ylabel "<Vi(0).Vj(t)>, (ang/ps)^2"
plot "FC_12" using 1:2 
pause -1 "Hit Return To Continue"

set title " FC_21"
set xlabel "time, ps"
set ylabel "<Vi(0).Vj(t)>, (ang/ps)^2"
plot "FC_21" using 1:2 
pause -1 "Hit Return To Continue"

set title " FC_22"
set xlabel "time, ps"
set ylabel "<Vi(0).Vj(t)>, (ang/ps)^2"
plot "FC_22" using 1:2 
pause -1 "Hit Return To Continue"

#set terminal postscript landscape color 18
#set output "temp.ps"
