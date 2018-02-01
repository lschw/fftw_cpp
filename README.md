# fftw_cpp
A c++ wrapper around the FFTW library for calculating Fourier transforms

Usage:

    // Create Fourier transform object with number of points 'N' and real space range 'tmax'
    FFT fft(N, tmax);
    
    // Perform Fourier transform on data array 'data' and store result in 'data_fft'
    fft.fft(data, data_fft);
    
    // Calculate corresponding frequencies and store in array 'f'
    fft.freq(f);
    
    // Shift order of frequencies and data to start with negative frequencies first
    fft.shift_freq(f, data_fft);


For a complete example see [test/test.cc](test/test.cc)
