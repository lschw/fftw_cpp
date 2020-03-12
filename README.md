# fftw_cpp
A c++ wrapper around the FFTW library for calculating 1d and 2d Fourier transforms

Usage:

    // Create Fourier transform object with number of points 'N' and real space range 'tmax'
    FFT fft(N, tmax);
    FFT2D fft2d(N1, N2, tmax1, tmax2);

    // Perform Fourier transform on data array 'data' and store result in 'data_fft'
    fft.fft(data, data_fft);
    fft2d.fft(data, data_fft);

    // Calculate corresponding frequencies and store in array 'f'
    fft.freq(f);
    fft2d.freq1(f1);
    fft2d.freq1(f2);

    // Shift order of frequencies and data to start with negative frequencies first
    fft.shift_freq(f, data_fft);
    fft2d.shift_freq(f1, f2, data_fft);


For a complete example see [test/test.cc](test/test.cc) or [test/test2d.cc](test/test2d.cc)
