#ifndef __FFTW_CPP__HH__
#define __FFTW_CPP__HH__

#include <cstring>
#include <vector>
#include <complex>
#include <fftw3.h>

typedef std::complex<double> dcomplex;
typedef std::vector<double> dvector;
typedef std::vector<dcomplex> dcvector;

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * Class representing a Fourier transform
 */
class FFT
{
    public:
        const size_t N; // Number of data points
        const double length; // Length of interval in real space
        const double sample_rate; // Sample rate (N/length)
        const double df; // (Angular) frequency step (2*pi/length)
    
    private:
        fftw_plan plan_fw;
        fftw_plan plan_bw;
    
    public:
        
        /**
         * Setup Fourier transform
         * @param N       Number of datapoints
         * @param length  Length of interval in real space
         */
        FFT(size_t N, double length) : N(N), length(length),
            sample_rate(N/length), df(2*M_PI/length)
        {
            #ifdef _OPENMP
            // Initialisize multithreaded FFT automatically if OpenMP is
            // available
            FFT::init_multithread(omp_get_max_threads());
            #endif
            
            plan_fw = fftw_plan_dft_1d(N, 0, 0, FFTW_FORWARD, FFTW_ESTIMATE);
            plan_bw = fftw_plan_dft_1d(N, 0, 0, FFTW_BACKWARD, FFTW_ESTIMATE);
        }
        
        
        /**
         * Clean up
         */
        ~FFT()
        {
            fftw_destroy_plan(plan_fw);
            fftw_destroy_plan(plan_bw);
        }
        
        
        /**
         * Calculate Fourier transform
         * @param in   Input data
         * @param out  Fourier transformed output data
         *             If in == out, the transformation is done in-place
         */
        void fft(dcvector& in, dcvector& out)
        {
            // Ensure in-place transformation
            if(in.data() != out.data()) {
                memcpy(out.data(), in.data(), N*sizeof(dcomplex));
            }
            fftw_execute_dft(plan_fw,
                reinterpret_cast<fftw_complex*>(out.data()),
                reinterpret_cast<fftw_complex*>(out.data())
            );
            
            // Scale amplitude as FFTW calculates unscaled coefficients
            for(size_t i = 0; i < N; ++i) {
                out[i] /= N;
            }
        }
        
        
        /**
         * Calculate inverse Fourier transform
         * @param in   Input data
         * @param out  Fourier transformed output data
         *             If in == out, the transformation is done in-place
         */
        void ifft(dcvector& in, dcvector& out)
        {
            // Ensure in-place transformation
            if(in.data() != out.data()) {
                memcpy(out.data(), in.data(), N*sizeof(dcomplex));
            }
            fftw_execute_dft(plan_bw,
                reinterpret_cast<fftw_complex*>(out.data()),
                reinterpret_cast<fftw_complex*>(out.data())
            );
        }
        
        
        /**
         * Calculate sample frequencies (angular frequency)
         * @param f  This array will store the frequency data. Format:
         *           [0, df, 2*df, ..., N/2*df,
         *            -(N/2-1)*df, -(N/2-2)*df, ..., -df]
         */
        void freq(dvector& f)
        {
            f.resize(N);
            for(size_t i = 0; i < N; ++i) {
                if(i <= N/2) {
                    // Positive frequencies first
                    f[i] = 2*M_PI*i*sample_rate/N;
                } else {
                    f[i] = -2*M_PI*(N-i)*sample_rate/N;
                }
            }
        }
        
        
        /**
         * Shift frequency and data array to order frequencies from negative
         * to positive
         * @param f     Frequency array
         * @param data  Data array
         */
        void shift_freq(dvector& f, dcvector& data)
        {
            dvector buf1(N);
            dcvector buf2(N);
            if(N%2 == 0) { // Even number of data points
                for(size_t i = 0; i < N/2+1; ++i) {
                    buf1[N/2-1+i] = f[i];
                    buf2[N/2-1+i] = data[i];
                    if(i < N/2-1) {
                        buf1[i] = f[N/2+1+i];
                        buf2[i] = data[N/2+1+i];
                    }
                }
            } else { // Odd number of data points
                buf1[N/2] = f[0];
                buf2[N/2] = data[0];
                for(size_t i = 0; i < N/2; ++i) {
                    buf1[N/2+1+i] = f[i+1];
                    buf2[N/2+1+i] = data[i+1];
                    buf1[i] = f[N/2+1+i];
                    buf2[i] = data[N/2+1+i];
                }
            }
            memcpy(f.data(), buf1.data(), N*sizeof(double));
            memcpy(data.data(), buf2.data(), N*sizeof(dcomplex));
        }
        
        
        /**
         * Initialisize FFTW to use multiple threads
         * Call this function before creating a FFT object
         * @param threads  Number of threads to use
         */
        static void init_multithread(int threads=4)
        {
            fftw_init_threads();
            fftw_plan_with_nthreads(threads);
        }
        
        
        /**
         * Uninitialisize multithreaded FFTW
         */
        static void clean_multithread()
        {
            fftw_cleanup_threads();
        }
};

#endif
