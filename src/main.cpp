// Author: Joshua Jansen van Vueren
// Date: 20/01/2021

#include <iostream>
#include "fft.h"
#include "stft.h"
#include "filter.h"
#include "doppler.h"
#include <cmath>
#include <complex>
#include <random>

using namespace std;

int main() {
    // stft vars
    int windowLen = 1024;
    int samplingFreq = 100000;
    int fftLen = 256;
    int signalLen = 50;
    int signalFreqStart = 10000;
    int signalFreq = 8000;

    // filter vars
    int filLen = 31;
    double atten = 10;
    double fLow = 6000;
    double fHigh = 13000;

    Filter fil; // create filter
    FFT fft;
    STFT stft(windowLen, samplingFreq, fftLen, false, "hamm");

    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution(0.0,0.9);

    // create a dummy sine wave with noise

    std::vector<std::complex<float>> signal(fftLen*signalLen,0);
    float m = (signalFreqStart - signalFreq)/(0.001);

    for(int i=0;i<fftLen*signalLen;i++) {
        float t = i * 1.0 / samplingFreq;
        if (t <= 0.1) {
            signal[i] = sinf(2* M_PI * (m*(t-0.1)*(t-0.1)+signalFreq) * t) + distribution(generator);
        } else {
            signal[i] = sinf(2* M_PI * signalFreq * t) + distribution(generator);
        } // else
        //cout << i << " @ " << t << " - " << signal[i] << endl;
    } // for

    std::vector<std::complex<float>> filter = fil.complexKaiserBesselFilterCoefficients(filLen, atten, fHigh, fLow, samplingFreq);
    std::vector<std::complex<float>> conv = fil.applyFilterByConv(signal,filter, signal.size(), filter.size());
    //std::vector<std::complex<float>> conv = signal;

    // truncates signal to multiple of fftLen
    stft.fitSignal(&conv);
    
    cout << "Full STFT" << endl;
    cout << ("---------------") << endl;
    stft.computeSTFT(&conv);
    
    //std::vector<std::vector<std::complex<float>>> stftResult = stft.getResult();
    std::vector<std::vector<float>> stftMagResult = stft.getMagResult();
    std::vector<float> freqBins = stft.getFreqBins();
    std::vector<float> timeBins = stft.getTimeBins();

    //for(int i=0;i<1;i++) { //timeBins.size()
    //    for (int j = 0; j < freqBins.size(); j++) {
    //        cout << timeBins[i] << " - " << freqBins[j] << " - " << stftMagResult[i][j] << endl;
    //    } // for
    //} // for

    // extract the principal component from the STFT
    std::vector<float> principle = fil.getPrinciple(stftMagResult,freqBins);

    for (int j = 0; j < principle.size(); j++) {
        cout << timeBins[j] << " - " << principle[j] << endl;
    } // for

    // TODO:
    // estimate doppler shift
    int transmitFreq = 10000;
    Doppler dop(transmitFreq);

    std::vector<float> projVel = dop.measureProjectileVelocity(principle);
    std::vector<float> projDis = dop.estimDistanceTravelled(projVel, timeBins);
    std::vector<float> transDis = dop.measureTransverseDistance(projVel, timeBins);

    cout << "-----------------------" << endl;
    cout << "------- velocity ------" << endl;
    cout << "-----------------------" << endl;
    for (int j = 0; j < projVel.size(); j++) {
        cout << timeBins[j] << " - " << projVel[j] << endl;
    } // for

    cout << "-----------------------" << endl;
    cout << "------- distance ------" << endl;
    cout << "-----------------------" << endl;
    for (int j = 0; j < projVel.size(); j++) {
        cout << timeBins[j] << " - " << projDis[j] << endl;
    } // for

    cout << "-----------------------" << endl;
    cout << "------ transverse -----" << endl;
    cout << "-----------------------" << endl;
    for (int j = 0; j < projVel.size(); j++) {
        cout << timeBins[j] << " - " << transDis[j] << endl;
    } // for
}

/***************************************************************
 *  OLD TEMP CODE

// create a dummy sine wave
    int signalFreq = 1000;
    float sampleFreq = 8000;
    int len = 8;
    float points[len];

    for(int i=0;i<len;i++) {
        float t = i * 1.0 / sampleFreq;
        points[i] = sinf(2* M_PI * signalFreq * t);
        cout << i << " @ " << t << " - " << points[i] << endl;
    }

    // create a copy of points for FFT calculation
    float result[len];
    std::copy ( &points[0], &points[len], &result[0] );
    FFT fft;

    std::complex<float> *radix_2_fft;
    radix_2_fft = fft.computeDitFft(&result[0],len);
    
    for(int i=0;i<len;i++) {
        cout << i << " - " << *radix_2_fft << endl;
        radix_2_fft++;
    }

    // convolution
    double seq1[] = {5,4,3,2,1};
    double seq2[] = {1,1,1,1};
    double * conv;
    conv = fil.ApplyFilterByConv(&seq1[0],&seq2[0], 5, 4);

    //for (int i = 0; i < 9; i++) {
    //    cout << i << " " << *(conv + i) << endl;
    //}
********************************************************/