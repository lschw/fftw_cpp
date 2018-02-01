GCC     = g++ -std=c++11
LIBS    = -lm -lpthread -lfftw3 -lfftw3_threads -fopenmp
CFLAGS  = -Wall -Wunused -Wextra -pedantic -Ofast

test: test.cc
	@$(GCC) -o test $(CFLAGS) test.cc $(LIBS)

run:
	./test
	gnuplot plot_data.gp
	gnuplot plot_data_fft.gp

clean:
	rm -f test data.dat data_fft.dat data.pdf data_fft.pdf
