#ifndef FILTER_H
#define FILTER_H

#include <complex>
#include <map>

class Filter {
    public:
        Filter();
        ~Filter();
        
        /**
         * @brief Estimate of the modified zero order Bessel function
         * sourced from https://www.atnf.csiro.au/computing/software/gipsy/
         * 
         * @param x Position at which to retrieve the Bessel estimate
         * @return double 
         */
        double bessi0(double x);

        /**
         * @brief Estimate of the modified zero order Bessel function
         * sourced from https://www.arc.id.au/FilterDesign.html
         * 
         * @param x Position at which to retrieve the Bessel estimate
         * @return double 
         */
        double bess_2(double x);

        /**
         * @brief Given an attenuation, this function returns the window shape
         * factor, which determines the trade off between main-lobe width and side lobe
         * level. See https://en.wikipedia.org/wiki/Kaiser_window
         * 
         * @param att Attenuation (dB) in side lobes
         * @return double 
         */
        double kaiserBesselWindowShape(int att);

        /**
         * @brief Returns an odd symmetrical modified sinc function of length len.
         * Where fHIgh and fLow are cutoff frequencies
         * 
         * @param fHigh High frequency cutoff (Hz)
         * @param fLow Low frequency cutoff (Hz)
         * @param fs Sampling frequency (Hz)
         * @param len Length of sinc function
         * @return double* 
         */
        double *sincFunction(double fHigh, double fLow, int fs, int len);

        /**
         * @brief Returns the time domain filter coefficients for a Kaiser-Bessel
         * bandpass filter - adapted from: https://www.arc.id.au/FilterDesign.html
         * 
         * @param len Length of filter
         * @param att Attenuation (dB) in side lobes
         * @param fHigh High frequency cutoff (Hz)
         * @param fLow Low frequency cutoff (Hz)
         * @param fs Sampling frequency (Hz)
         * @return double* 
         */
        double *kaiserBesselFilterCoefficients(int len,double att, double fHigh, double fLow, int fs);

        /**
         * @brief Returns the time domain filter coefficients for a Kaiser-Bessel
         * bandpass filter - adapted from: https://www.arc.id.au/FilterDesign.html
         * 
         * @param len Length of filter
         * @param att Attenuation (dB) in side lobes
         * @param fHigh High frequency cutoff (Hz)
         * @param fLow Low frequency cutoff (Hz)
         * @param fs Sampling frequency (Hz)
         * @return std::vector<std::complex<float>> 
         */
        std::vector<std::complex<float>> complexKaiserBesselFilterCoefficients(int len, double att, double fHigh, double fLow, int fs);

        /**
         * @brief Returns an estimated n point hamming
         * 
         * @param n 
         * @return std::vector<float> 
         */
        std::vector<float> hammingWindow(int n);

        /**
         * @brief Returns an estimated n point hamming
         * 
         * @param n 
         * @return std::vector<std::complex<float>> 
         */
        std::vector<std::complex<float>> complexHammingWindow(int n);

        /**
         * @brief Convolves two signals to produce another of length len1+ len2
         * 
         * @param seq First sequence to convolve
         * @param fil Second sequence to convolve
         * @param seqLen First sequence length
         * @param filLen Second sequence length
         * @return double* 
         */
        double *applyFilterByConv(double * seq, double * fil, int seqLen, int filLen);

        /**
         * @brief Convolves two signals to produce another of length len1+ len2
         * 
         * @param seq First sequence to convolve
         * @param fil Second sequence to convolve
         * @param seqLen First sequence length
         * @param filLen Second sequence length
         * @return std::vector<std::complex<float>> 
         */
        std::vector<std::complex<float>> applyFilterByConv(std::vector<std::complex<float>> seq, std::vector<std::complex<float>> fil, int seqLen, int filLen);

        /**
         * @brief Get the principle frequency (argmax) given a 2d STFT and the
         * frequency bins
         * 
         * @param inp 2D STFT
         * @param arg Frequency Bins
         * @return std::vector<float> 
         */
        std::vector<float> getPrinciple(std::vector<std::vector<float>> inp, std::vector<float> arg);
};

#endif // FILTER_H