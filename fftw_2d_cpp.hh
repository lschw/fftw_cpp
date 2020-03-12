#ifndef __FFTW2D_CPP__HH__
#define __FFTW2D_CPP__HH__

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
 * Class representing a 2D Fourier transform
 */
class FFT2D
{
    public:
        const size_t N1; // Number of data points 1
        const double length1; // Length of interval in real space 1
        const double sample_rate1; // Sample rate (N/length) 1
        const double df1; // (Angular) frequency step (2*pi/length) 1
        const size_t N2; // Number of data points 2
        const double length2; // Length of interval in real space 2
        const double sample_rate2; // Sample rate (N/length) 2
        const double df2; // (Angular) frequency step (2*pi/length) 2
        const size_t N; // Total number of data points = N1*N2

    private:
        fftw_plan plan_fw;
        fftw_plan plan_bw;

    public:

        /**
         * Setup Fourier transform
         * @param N1       Number of datapoints in first dim.
         * @param N2       Number of datapoints in second dim.
         * @param length1  Length of interval in real space in first dim.
         * @param length2  Length of interval in real space in second dim.
         */
        FFT2D(size_t N1, size_t N2, double length1, double length2) :
            N1(N1), length1(length1),
            sample_rate1(N1/length1), df1(2*M_PI/length1),
            N2(N2), length2(length2),
            sample_rate2(N2/length2), df2(2*M_PI/length2),
            N(N1*N2)
        {
            #ifdef _OPENMP
            // Initialisize multithreaded FFT automatically if OpenMP is
            // available
            FFT2D::init_multithread(omp_get_max_threads());
            #endif

            plan_fw = fftw_plan_dft_2d(
                N1, N2, 0, 0, FFTW_FORWARD, FFTW_ESTIMATE);
            plan_bw = fftw_plan_dft_2d(
                N1, N2, 0, 0, FFTW_BACKWARD, FFTW_ESTIMATE);
        }


        /**
         * Clean up
         */
        ~FFT2D()
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
        void _freq(dvector& f, size_t N, double sr)
        {
            f.resize(N);
            for(size_t i = 0; i < N; ++i) {
                if(i <= N/2) {
                    // Positive frequencies first
                    f[i] = 2*M_PI*i*sr/N;
                } else {
                    f[i] = -2*M_PI*(N-i)*sr/N;
                }
            }
        }
        void freq1(dvector& f)
        {
            return _freq(f, N1, sample_rate1);
        }
        void freq2(dvector& f)
        {
            return _freq(f, N2, sample_rate2);
        }


        /**
         * Shift frequency and data array to order frequencies from negative
         * to positive
         * @param f1     Frequency array1
         * @param f2     Frequency array2
         * @param data  Data array
         */
        void shift_freq(dvector& f1, dvector& f2, dcvector& data)
        {
            dvector buf1(N1);
            dvector buf2(N2);
            dcvector bufd(N);

            // Shift first dimension
            if(N1%2 == 0) { // Even number of data points
                for(size_t i = 0; i < N1/2+1; ++i) {
                    buf1[N1/2-1+i] = f1[i];

                    for(size_t j = 0; j < N2; ++j) {
                        bufd[(N1/2-1+i)*N2 + j] = data[i*N2 + j];
                    }
                    if(i < N1/2-1) {
                        buf1[i] = f1[N1/2+1+i];
                        for(size_t j = 0; j < N2; ++j) {
                            bufd[i*N2 + j] = data[(N1/2+1+i)*N2 + j];
                        }
                    }
                }
            } else { // Odd number of data points
                buf1[N1/2] = f1[0];
                bufd[N1/2] = data[0];
                for(size_t i = 0; i < N1/2; ++i) {
                    buf1[N1/2+1+i] = f1[i+1];

                    for(size_t j = 0; j < N2; ++j) {
                        bufd[(N1/2+1+i)*N2 + j] = data[(i+1)*N2 + j];
                    }
                    buf1[i] = f1[N1/2+1+i];
                    for(size_t j = 0; j < N2; ++j) {
                        bufd[i*N2 + j] = data[(N1/2+1+i)*N2 + j];
                    }
                }
            }

            // Shift second dimension
            if(N2%2 == 0) { // Even number of data points
                for(size_t i = 0; i < N2/2+1; ++i) {
                    buf2[N2/2-1+i] = f2[i];

                    for(size_t j = 0; j < N1; ++j) {
                        data[j*N2 + N2/2-1+i] = bufd[j*N2 + i];
                    }
                    if(i < N2/2-1) {
                        buf2[i] = f2[N2/2+1+i];
                        for(size_t j = 0; j < N1; ++j) {
                            data[j*N2 + i] = bufd[j*N2 + N2/2+1+i];
                        }
                    }
                }
            } else { // Odd number of data points
                buf2[N2/2] = f2[0];

                for(size_t j = 0; j < N1; ++j) {
                    data[j*N2 + N2/2] = bufd[j*N2 + 0];
                }
                for(size_t i = 0; i < N2/2; ++i) {
                    buf2[N2/2+1+i] = f2[i+1];

                    for(size_t j = 0; j < N1; ++j) {
                        data[j*N2 + N2/2+1+i] = bufd[j*N2 + i+1];
                    }
                    buf2[i] = f2[N2/2+1+i];
                    for(size_t j = 0; j < N1; ++j) {
                        data[j*N2 + i] = bufd[j*N2 + N2/2+1+i];
                    }
                }
            }

            memcpy(f1.data(), buf1.data(), N1*sizeof(double));
            memcpy(f2.data(), buf2.data(), N2*sizeof(double));
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
