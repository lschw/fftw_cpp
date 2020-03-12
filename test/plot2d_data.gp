set terminal pdf color size 4,3
set output "data2d.pdf"

set xlabel "t1"
set ylabel "t2"
set zlabel "f(t1,t2)"

set xrange[0:1]
set yrange[0:1]
set view map
splot "data2d.dat" u 1:2:3 w pm3d title ""
