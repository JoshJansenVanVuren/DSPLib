#include <iostream>
#include <complex>
#include <string>
#include "fft.h"
#include "stft.h"
#include "filter.h"

using namespace std;

STFT::STFT(int windowLen, int samplingFreq, int fftLen, bool ignoreNquist, std::string window) :
windowLen(windowLen),samplingFreq(samplingFreq),fftLen(fftLen),ignoreNquist(ignoreNquist),window(window) {
    FFT::setFFTLen(fftLen);

    if (windowLen > fftLen) {
        STFT::zeroPadding = true;
    } // if

    if (ignoreNquist) {
        for (int i = 0; i < windowLen; i++) {
            STFT::freqBins.push_back(i*(samplingFreq/windowLen));
        } // for
    } else {
        for (int i = 0; i < windowLen/2; i++) {
            STFT::freqBins.push_back(i*(samplingFreq/windowLen));
        } // for
    } // else
}

STFT::~STFT() { }

void STFT::fitSignal(std::vector<std::complex<float>> *signal) {
    int lastIndex = (*signal).size() - (((*signal).size() / STFT::fftLen) * STFT::fftLen);

    for (int i = 0; i < lastIndex; i++) {
        (*signal).pop_back();
    }
}

void STFT::computeSTFT(std::vector<std::complex<float>> *signal) {
    // the assumption is made here that the
    // signal has been cleaned and padded
    // i.e. length of input % windowLen == 0

    // calculate the time bins
    for (int i = 0; i < (*signal).size(); i+=fftLen) {
        STFT::timeBins.push_back(((float) i)/(STFT::samplingFreq));
    }

    // compute twiddles for all FFTs
    FFT::computeTwiddles();

    int len = FFT::getFFTLen();    

    for (int n = 0; n < (*signal).size(); n+=len) {
        vector<std::complex<float>> subSignal((*signal).begin() + n, (*signal).begin() + n + len);      

        // add zero padded tokens
        if (STFT::zeroPadding) {
            for (int i = 0; i < (windowLen - len); i++) {
                subSignal.push_back(0);
            } // for
        } // if

        // add the window function
        if (STFT::window == "hamm") {
            vector<std::complex<float>> window = Filter::complexHammingWindow(len+1);

            for (int i = 0; i < len; i++) {
                subSignal[i] *= window[i];
            } // for
        } // for
        

        // decimate signal in time
        FFT::computeDit(&subSignal,windowLen);

        // compute butterflies
        for (int n = 2; n <= windowLen; n*=2) {
            for (int k = 0; k < windowLen/n; k++) {
                STFT::nPointButterfly(&subSignal,k,n);
            } // for
        } // for

        if (ignoreNquist) {
            STFT::result.push_back(subSignal);
        } else {
            vector<std::complex<float>> halfSignal(subSignal.begin(), subSignal.begin() + windowLen/2);
            STFT::result.push_back(halfSignal);
        } // else
    } // for
}

std::vector<float> STFT::getFreqBins() {
    return STFT::freqBins;
}
        
std::vector<float> STFT::getTimeBins() {
    return STFT::timeBins;
}

std::vector<std::vector<std::complex<float>>> STFT::getResult() {
    return STFT::result;
}

std::vector<std::vector<float>> STFT::getMagResult() {
    std::vector<std::vector<float>> abs;
    for (int i = 0; i < STFT::result.size(); i++) {
        std::vector<float> abs_vec;
        for (int j = 0; j < STFT::result[i].size(); j++) {
            abs_vec.push_back(std::abs(STFT::result[i][j]));
        } // for
        abs.push_back(abs_vec);
    } // for
    
    return abs;
}