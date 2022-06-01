# data
if (!exists("m")) m="matriz-sistema.dat"
set terminal x11 0
set nokey

# labeling
# set xtics 5 out nomirror
# set ytics 5 out nomirror
set xlabel "x"
set ylabel "y"
set title m

# autoescalado
set autoscale xfix
set autoscale yfix
set autoscale cbfix

# paleta de color
set palette gray negative

# plot
plot m matrix with image

# save (png)
set term png             
set output m[1:strlen(m)-4].".png" 
replot
set term x11
