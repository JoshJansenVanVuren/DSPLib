#ifndef STFT_H
#define STFT_H

#include <complex>
#include <map>
#include <vector>
#include <string>
#include "fft.h"
#include "filter.h"

class STFT : public FFT, public Filter {
    public:
        /**
         * @brief Construct a new STFT object
         * 
         * @param windowLen Length of single FFT (if larger than fftLen, sequence is zero padded) Must be power of two.
         * @param samplingFreq Sampling frequency (Hz)
         * @param fftLen Length of FFT. Must be power of two.
         * @param ignoreNquist If Nquist is ignored the full spectrum is returned, otherwise spectrum is samplingFreq/2
         * @param window Type of window. "hamm" or "none".
         */
        STFT(int windowLen, int samplingFreq, int fftLen, bool ignoreNquist, std::string window="hamm");
        
        ~STFT();

        /**
         * @brief This function truncates a signal of arbitrary length
         * to a power of 2 for the FFT
         * 
         * @param signal 
         */
        void fitSignal(std::vector<std::complex<float>> *signal);

        /**
         * @brief Computes the STFT of a given signal
         * 
         * @param signal 
         */
        void computeSTFT(std::vector<std::complex<float>> *signal);

        /**
         * @brief Get the result from the STFT computation, returns
         * a 2d complex vector
         * 
         * @return std::vector<std::vector<std::complex<float>>> 
         */
        std::vector<std::vector<std::complex<float>>> getResult();

        /**
         * @brief Get the magnitude result from the STFT computation,
         * returns a 2d float vector
         * 
         * @return std::vector<std::vector<float>> 
         */
        std::vector<std::vector<float>> getMagResult();

        /**
         * @brief Get the freq bins given the STFT parameters
         * 
         * @return std::vector<float> 
         */
        std::vector<float> getFreqBins();

        /**
         * @brief Get the time bins given the STFT parameters
         * 
         * @return std::vector<float> 
         */
        std::vector<float> getTimeBins();
    
    private:
        int windowLen;
        int samplingFreq;
        int fftLen;
        bool zeroPadding;
        bool ignoreNquist;
        std::string window;
        
        std::vector<std::vector<std::complex<float>>> result;
        std::vector<float> freqBins;
        std::vector<float> timeBins;
};

#endif // STFT_H