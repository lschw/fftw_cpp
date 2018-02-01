#include <fstream>

#include "../fftw_cpp.hh"

int main()
{
    size_t N = 1000;
    dcvector data(N);
    dcvector data_fft(N);
    dvector t(N);
    dvector f(N);
    
    // Create some test data
    double w1 = 2*M_PI;
    double w2 = 2*M_PI*3;
    double tmax = 50;
    for(size_t i = 0; i < N; ++i) {
        t[i] = i*tmax/N;
        data[i] = sin(w1*t[i]) + sin(w2*t[i]);
    }
    
    // Fourier transform
    FFT fft(N, tmax);
    fft.fft(data, data_fft);
    fft.freq(f);
    fft.shift_freq(f, data_fft);
    
    // Save
    std::ofstream fh1;
    std::ofstream fh2;
    fh1.open("data.dat");
    fh2.open("data_fft.dat");
    fh1 << "# t Re[f(t)] Im[f(t)]\n";
    fh2 << "# f Re[f(w)] Im[f(w)]\n";
    for(size_t i = 0; i < N; ++i) {
        fh1 << t[i] << " ";
        fh1 << data[i].real() << " " << data[i].imag()<< "\n";
        fh2 << f[i] << " ";
        fh2 << data_fft[i].real() << " " << data_fft[i].imag()<< "\n";
    }
    fh1.close();
    fh2.close();
} 
