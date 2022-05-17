# data
if (!exists("m")) m="resultados/energia-especifica.dat"
set terminal x11 0
set nokey

# labeling
set xlabel "Iteraciones"
set ylabel "Energia especifica"
set title "Evolucion de energia especifica"

# autoescalado
set autoscale xfix
set autoscale yfix

# Set linestyles
set style line 1 \
    linecolor rgb '#0060ad' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.5
set style line 2 \
    linecolor rgb '#8df21b' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.5
set style line 3 \
    linecolor rgb '#de2e2e ' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.5

# plot
plot m with lines linestyle 3

# save (png)
set term png             
set output m[1:strlen(m)-4].".png" 
replot
set term x11