GCC     = g++ -std=c++11
LIBS    = -lm -lpthread -lfftw3 -lfftw3_threads -fopenmp
CFLAGS  = -Wall -Wunused -Wextra -pedantic -Ofast

all: test test2d

test: test.cc
	@$(GCC) -o test $(CFLAGS) test.cc $(LIBS)

test2d: test2d.cc
	@$(GCC) -o test2d $(CFLAGS) test2d.cc $(LIBS)

run:
	./test
	./test2d
	gnuplot plot_data.gp
	gnuplot plot_data_fft.gp
	gnuplot plot2d_data.gp
	gnuplot plot2d_data_fft.gp

clean:
	rm -f test test2d data.dat data_fft.dat data.pdf data_fft.pdf data2d.dat data2d_fft.dat data2d.pdf data2d_fft.pdf
