#include <iostream>
#include <complex>
#include <cmath>
#include <map>
#include "fft.h"

using namespace std;

FFT::FFT() {}
FFT::~FFT() {}

void FFT::setFFTLen(int l) {
    FFT::len = l;
}

int FFT::getFFTLen() {
    return FFT::len;
}

std::complex<float> *FFT::getResult() {
    return FFT::radix_2_fft;
}

void FFT::computeDit(float *points,int len) {
    // base case is len == 2
    if (len == 2) {
        return;
    } else {
        // place even index in first half
        // and odd indexes in second
        float * even_temp = (float*) malloc(len/2 * sizeof(float*));
        float * odd_temp = (float*) malloc(len/2 * sizeof(float*));

        for (int i = 0; i < len/2; i++) {
            float temp = *(points+i*2);
            *(even_temp+i) = temp;

            temp = *(points+1+i*2);
            *(odd_temp+i) = temp;
        } // for

        for (int i = 0; i < len/2; i++) {
            float temp = *(even_temp+i);
            *(points+i) = temp;

            temp = *(odd_temp+i);
            *(points+len/2+i) = temp;
        } // for
        
        // call two len / 2 decimation functions
        FFT::computeDit(points,len/2);
        FFT::computeDit(points+len/2,len/2);

    } // else
}

void FFT::computeDit(std::vector<std::complex<float>> *signal, int len) {
    for (int i = len; i > 2; i/=2) {
        std::vector<std::complex<float>> temp (len,0);
        for (int j = 0; j < len; j+=i) {
            // place even index in first half
            // and odd indexes in second
            for (int k = 0; k < i/2; k++) {
                temp[k + j] = (*signal)[k*2 + j];
                temp[k + j + i/2] = (*signal)[k*2 + j + 1];
            } // for
        } // for

        (*signal).swap(temp);
    } // for
}

void FFT::computeTwiddles() {
    // we want to instantiate all the possibilities for W_N^k
    // std::map<std::pair<int,int>,std::complex<float>> twiddles;
    FFT::twiddles.clear();

    for (int k = 0; k < FFT::len/2; k++) {
        complex<float> twid = std::polar(((float) 1), ((float)(-(2*M_PI*k)/len)));
        
        FFT::twiddles[k] = twid;
    } // for
}

std::complex<float> FFT::getTwiddle(int k) {
    return FFT::twiddles[k];
}

void FFT::toComplex(float *points,int len) {
    // allocate memory for fft result
    FFT::radix_2_fft = (std::complex<float>*) malloc(len * sizeof(std::complex<float>*));

    // cast float -> complex
    for (int i = 0; i < len; i++) {
        complex<float> tmp (*(points+i),(float)0);
        *(FFT::radix_2_fft + i) = tmp;
    } // for
}

std::complex<float> FFT::complexMul(std::complex<float> val1, std::complex<float> val2) {
    std::complex<float> result;

    result.real((val1.real() * val2.real()) - (val1.imag() * val2.imag()));
    result.imag((val1.real() * val2.imag()) + (val1.imag() * val2.real()));

    return result;
}

void FFT::nPointButterfly(std::complex<float> *points, int n) {
    // add twiddle factors to second half of signal
    for (int i = 0; i < n/2; i++) {
        std::complex<float> mul;
        std::complex<float> twid = FFT::getTwiddle(i*FFT::len/n);
        mul = FFT::complexMul(*(points + n/2 + i),twid);
        *(points + n/2 + i) = mul;
    } // for

    // apply cross over multiplications
    for (int i = 0; i < n/2; i++) {
        std::complex<float> first_temp;
        std::complex<float> second_temp;

        first_temp =  *(points + i) + *(points + i + n/2);
        second_temp = *(points + i) - *(points + i + n/2);
        
        *(points + i) = first_temp;
        *(points + i + n/2) = second_temp;
    } // for
}

void FFT::nPointButterfly(std::vector<std::complex<float>> *signal, int k, int n) {
    // add twiddle factors to second half of signal
    for (int i = 0; i < n/2; i++) {
        std::complex<float> twid = FFT::getTwiddle(i*FFT::len/n);
        (*signal)[i + k*n + n/2] *= twid;
    } // for

    // apply cross over computation
    for (int i = 0; i < n/2; i++) {
        std::complex<float> first_temp;
        std::complex<float> second_temp;

        first_temp =  (*signal)[i + k*n] + (*signal)[i + k*n + n/2];
        second_temp = (*signal)[i + k*n] - (*signal)[i + k*n + n/2];
        
        (*signal)[i + k*n] = first_temp;
        (*signal)[i + k*n + n/2] = second_temp;
    } // for
}

std::complex<float> *FFT::computeDitFft(float *points,int len) {
    FFT::len = len;

    // precompute the twiddles for the DFT
    FFT::computeTwiddles();

    // decimate the signal in time
    FFT::computeDit(points,len);

    // create complex input
    FFT::toComplex(points,len);

    // reconstruct the signal in z-domain
    for (int n = 2; n <= len; n*=2) {
        for (int k = 0; k < len/n; k++) {
            FFT::nPointButterfly((FFT::radix_2_fft+k*n),n);
        } // for
    } // for
    clog << "fft complete" << endl;

    return FFT::radix_2_fft;
}