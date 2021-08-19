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

set xlabel "CF4 (2) molec/uc"
set ylabel "dlnf1/dlnc1"
#set key left
set xrange [0:48]
set title " dlnf1/dlnc1, Methane=1, CF4=2"
plot "Me8.dat" using 2:7 title "Methane =8/uc", "Me16.dat" using 2:7 title "Methane =16/uc", "Me24.dat" using 2:7 title "Methane =24/uc", "Me32.dat" using 2:7 title "Methane =32/uc", "Me40.dat" using 2:7 title "Methane =40/uc" 
pause -1 "Hit Return To Continue"


set xlabel "CF4 (2) molec/uc"
set ylabel "dlnf1/dlnc2"
#set key left
set xrange [0:48]
set title " dlnf1/dlnc2, Methane=1, CF4=2"
plot "Me8.dat" using 2:8 title "Methane =8/uc", "Me16.dat" using 2:8 title "Methane =16/uc", "Me24.dat" using 2:8 title "Methane =24/uc", "Me32.dat" using 2:8 title "Methane =32/uc", "Me40.dat" using 2:8 title "Methane =40/uc" 
pause -1 "Hit Return To Continue"


set xlabel "CF4 (2) molec/uc"
set ylabel "dlnf2/dlnc2"
#set key left
set xrange [0:48]
set title " dlnf2/dlnc2, Methane=1, CF4=2"
plot "Me8.dat" using 2:9 title "Methane =8/uc", "Me16.dat" using 2:9 title "Methane =16/uc", "Me24.dat" using 2:9 title "Methane =24/uc", "Me32.dat" using 2:9 title "Methane =32/uc", "Me40.dat" using 2:9 title "Methane =40/uc" 
pause -1 "Hit Return To Continue"

set xlabel "CF4 (2) molec/uc"
set ylabel "dlnf2/dlnc1"
#set key left
set xrange [0:48]
set title " dlnf2/dlnc1, Methane=1, CF4=2"
plot "Me8.dat" using 2:10 title "Methane =8/uc", "Me16.dat" using 2:10 title "Methane =16/uc", "Me24.dat" using 2:10 title "Methane =24/uc", "Me32.dat" using 2:10 title "Methane =32/uc", "Me40.dat" using 2:10 title "Methane =40/uc" 
pause -1 "Hit Return To Continue"

set terminal postscript landscape color 18

set output "f1c1.eps"
set xlabel "CF4 (2) molec/uc"
set ylabel "dlnf1/dlnc1"
set xrange [0:48]
set title " dlnf1/dlnc1, Methane=1, CF4=2"
plot "Me8.dat" using 2:7 title "Methane =8/uc", "Me16.dat" using 2:7 title "Methane =16/uc", "Me24.dat" using 2:7 title "Methane =24/uc", "Me32.dat" using 2:7 title "Methane =32/uc", "Me40.dat" using 2:7 title "Methane =40/uc" 

set output "f1c2.eps"
set xlabel "CF4 (2) molec/uc"
set ylabel "dlnf1/dlnc2"
set xrange [0:48]
set title " dlnf1/dlnc2, Methane=1, CF4=2"
plot "Me8.dat" using 2:8 title "Methane =8/uc", "Me16.dat" using 2:8 title "Methane =16/uc", "Me24.dat" using 2:8 title "Methane =24/uc", "Me32.dat" using 2:8 title "Methane =32/uc", "Me40.dat" using 2:8 title "Methane =40/uc" 

set output "f2c2.eps"
set xlabel "CF4 (2) molec/uc"
set ylabel "dlnf2/dlnc2"
set xrange [0:48]
set title " dlnf2/dlnc2, Methane=1, CF4=2"
plot "Me8.dat" using 2:9 title "Methane =8/uc", "Me16.dat" using 2:9 title "Methane =16/uc", "Me24.dat" using 2:9 title "Methane =24/uc", "Me32.dat" using 2:9 title "Methane =32/uc", "Me40.dat" using 2:9 title "Methane =40/uc" 

set output "f2c1.eps"
set xlabel "CF4 (2) molec/uc"
set ylabel "dlnf2/dlnc1"
set xrange [0:48]
set title " dlnf2/dlnc1, Methane=1, CF4=2"
plot "Me8.dat" using 2:10 title "Methane =8/uc", "Me16.dat" using 2:10 title "Methane =16/uc", "Me24.dat" using 2:10 title "Methane =24/uc", "Me32.dat" using 2:10 title "Methane =32/uc", "Me40.dat" using 2:10 title "Methane =40/uc" 

