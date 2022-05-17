# data
if (!exists("m")) m="resultados/masa.dat"
set terminal x11 0
set nokey

# labeling
set xlabel "Iteraciones"
set ylabel "Masa"
set title "Evolucion de la masa"

# autoescalado
set autoscale xfix
set autoscale yfix

# Set linestyles
set style line 1 \
    linecolor rgb '#ffc300' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.5

# plot
plot m with lines linestyle 1

# save (png)
set term png             
set output m[1:strlen(m)-4].".png" 
replot
set term x11