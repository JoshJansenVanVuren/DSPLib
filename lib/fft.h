#ifndef FFT_H
#define FFT_H

#include <complex>
#include <vector>
#include <map>

class FFT {
    public:
        FFT();
        ~FFT();

        /**
         * @brief Inplace recursive function to decimate a sequence for radix-2 FFT
         * 
         * @param points Pointer to sequence to decimate
         * @param len The length of the array
         * @return
         */
        void computeDit(float *points, int len);

        /**
         * @brief Inplace recursive function to decimate a sequence for radix-2 FFT.
         * 
         * @param signal Signal to decimate
         * @param len Length of decimation
         */
        void computeDit(std::vector<std::complex<float>> *signal, int len);
        
        /**
         * @brief Set the FFT objects length
         * 
         * @param len 
         */
        void setFFTLen(int len);

        /**
         * @brief Gets the FFT objects length
         * 
         * @return int 
         */
        int getFFTLen();

        /**
         * @brief Get the computed FFT
         * 
         * @return std::complex<float>* 
         */
        std::complex<float> *getResult();
       
        /**
         * @brief  Initialises an FFT objects twiddles
         * 
         */
        void computeTwiddles();

        /**
         * @brief Get a value from the twiddle map
         * 
         * @param k 
         * 
         * @return std::complex<float> 
         */
        std::complex<float> getTwiddle(int k);

        /**
         * @brief Casts values from a float pointer to a complex pointer (onto member variable radix_2_fft)
         * 
         * @param points 
         * @param len 
         */
        void toComplex(float *points,int len);

        /**
         * @brief Computes the complex mutliplication of two signals by 
         * decomposing into 4 real multiplications
         * 
         * @param val1 
         * @param val2 
         * @return std::complex<float> 
         */
        std::complex<float> complexMul(std::complex<float> val1, std::complex<float> val2);

        /**
         * @brief Computes an n point butterfly multiplication on signal
         * 
         * @param points Signal to compute the butterfly on
         * @param len Length of the butterfly
         */
        void nPointButterfly(std::complex<float> *points, int len);

        /**
         * @brief Computes an n point butterfly multiplication on signal
         * 
         * @param signal Signal to compute the butterfly on
         * @param k Start index for the butterfly
         * @param n Length of the butterfly
         */
        void nPointButterfly(std::vector<std::complex<float>> *signal, int k, int n);

        /**
         * @brief Computes the Cooley–Tukey FFT algorithm FFT of a real DIT sequence
         * 
         * @param points Signal from which the FFT is computed
         * @return float* 
         */
        std::complex<float> *computeDitFft(float *points, int len);
    
    private:
        int len;
        std::complex<float> *radix_2_fft;
        std::map<int,std::complex<float>> twiddles;
};

#endif // FFT_H