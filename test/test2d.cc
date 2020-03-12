#include <iostream>
#include <fstream>
#include "../fftw_2d_cpp.hh"

int main()
{
    size_t N1 = 100;
    size_t N2 = 500;
    size_t N = N1*N2;
    dcvector data(N);
    dcvector data_fft(N);
    dvector t1(N1);
    dvector t2(N2);
    dvector f1(N1);
    dvector f2(N2);

    // Create test data
    double w1 = 2*M_PI;
    double w2 = 2*M_PI*3;
    double w3 = 2*M_PI*5;
    double w4 = 2*M_PI*7;
    double tmax1 = 5;
    double tmax2 = 10;
    for(size_t i = 0; i < N1; ++i) {
        t1[i] = i*tmax1/N1;
        for(size_t j = 0; j < N2; ++j) {
            t2[j] = j*tmax2/N2;
            data[i*N2 + j] = (sin(w1*t1[i]) + sin(w2*t1[i]))
                * (sin(w3*t2[j]) + sin(w4*t2[j]));
        }
    }

    // Fourier transform
    FFT2D fft(N1, N2, tmax1, tmax2);
    fft.fft(data, data_fft);
    fft.freq1(f1);
    fft.freq2(f2);
    fft.shift_freq(f1, f2, data_fft);

    // Save
    std::ofstream fh1;
    std::ofstream fh2;
    fh1.open("data2d.dat");
    fh2.open("data2d_fft.dat");
    fh1 << "# t1 t2 Re[f(t1,t2)] Im[f(t1,t2)]\n";
    fh2 << "# f1 f2 Re[f(f1,f2)] Im[f(f1,f2)]\n";
    for(size_t i = 0; i < N1; ++i) {
        for(size_t j = 0; j < N2; ++j) {
            fh1 << t1[i] << " ";
            fh1 << t2[j] << " ";
            fh1 << data[i*N2+j].real() << " ";
            fh1 << data[i*N2+j].imag() << "\n";
            fh2 << f1[i] << " ";
            fh2 << f2[j] << " ";
            fh2 << data_fft[i*N2+j].real() << " ";
            fh2 << data_fft[i*N2+j].imag() << "\n";
        }
        fh1 << "\n";
        fh2 << "\n";
    }
    fh1.close();
    fh2.close();
}
