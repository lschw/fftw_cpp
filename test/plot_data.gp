set terminal pdf color size 4,3
set output "data.pdf"

set xlabel "t"
set ylabel "f(t)"

set xrange[0:10]
plot "data.dat" u 1:2 w l title ""
