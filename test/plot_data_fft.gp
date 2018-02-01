set terminal pdf color size 4,3
set output "data_fft.pdf"

set xlabel "f"
set ylabel "Power spectrum"

set xrange[0:5]
plot "data_fft.dat" u ($1/(2*pi)):($2*$2+$3*$3) w l title ""
