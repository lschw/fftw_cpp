set terminal pdf color size 4,3
set output "data2d_fft.pdf"

set xlabel "f1"
set ylabel "f2"
set zlabel "Power spectrum"

set xrange[0:5]
set yrange[0:10]
set view map
splot "data2d_fft.dat" u ($1/(2*pi)):($2/(2*pi)):(($3*$3+$4*$4)*10) w pm3d title ""


