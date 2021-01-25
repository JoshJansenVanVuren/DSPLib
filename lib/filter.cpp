#include <iostream>
#include <cmath>
#include <vector>
#include "filter.h"

using namespace std;

Filter::Filter() {}
Filter::~Filter() {}

double Filter::bess_2(double x) {
  /*
   * This function calculates the zeroth order Bessel function
   */
  double d = 0, ds = 1, s = 1;
  do {
    d += 2;
    ds *= x*x/(d*d);
    s += ds;
  } while (ds > s*1e-6);
  return s;
}

double Filter::bessi0( double x ) {
    /**
     * @brief Sourced from https://www.atnf.csiro.au/computing/software/gipsy/
     * 
    /* PURPOSE: Evaluate modified Bessel function In(x) and n=0.  */
    /*------------------------------------------------------------*/
   double ax,ans;
   double y;

   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
   } else {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
         +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
         +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
         +y*0.392377e-2))))))));
   }
   return ans;
}

double Filter::kaiserBesselWindowShape(int att) {
    if (att > 50) {
        return 0.1102 * (((double) att) - 8.7);
    } else if ((att >= 21) && (att <= 50)) {
        return 0.5842 * pow((((double) att) - 21),0.4) + 0.07886 * (((double) att) - 21);
    } else {
        return 0;
    }
}

double* Filter::sincFunction(double fHigh, double fLow, int fs, int len) {
    double* result = (double*) malloc((len+1)/2 * sizeof(double*));

    double fHighFw = fHigh / ((double) fs);
    double fLowFw = fLow / ((double) fs);

    *result = ( 2 * (fHigh - fLow)) / (fs);
    
    for (int i = 1; i <= (len+1)/2; i++) {
        double temp = ( sin(2*M_PI*i*fHighFw) - sin(2*M_PI*i*fLowFw)) / (M_PI * i);
        *(result + i) = temp;
    } // for
    

    return result;
}

double * Filter::kaiserBesselFilterCoefficients(int len, double att, double fHigh, double fLow, int fs) {
    /**
     * Read more: https://en.wikipedia.org/wiki/Kaiser_window
     * 
     */

    double * result = (double*) malloc(len * sizeof(double*));
    double * sinc = Filter::sincFunction(fHigh, fLow, fs, len);
    int half_len = (len-1)/2;
    double alpha = Filter::kaiserBesselWindowShape(att);
    //double besselAlpha =  Filter::bessi0(alpha);
    double besselAlpha2 = Filter::bess_2(alpha);

    for (int i = 0; i <= half_len; i++) {
        double bessArg = alpha * sqrt(1 - ((double)i *i)/((double) half_len * half_len));
        //double besselI = Filter::bessi0(bessArg);
        double besselI2 = Filter::bess_2(bessArg);
        *(result + half_len + i) = (*(sinc + i) * besselI2) / besselAlpha2;
    } // for

    for (int i = 0; i < half_len; i++) {
        double temp = *(result + len - 1 - i);
        *(result + i) = temp;
    } // for

    return result;
}

std::vector<std::complex<float>> Filter::complexKaiserBesselFilterCoefficients(int len, double att, double fHigh, double fLow, int fs) {
    double * result = Filter::kaiserBesselFilterCoefficients(len, att, fHigh, fLow, fs);
    std::vector<std::complex<float>> complexResult;

    for (int i = 0; i < len; i++) {
        complexResult.push_back((float) *(result + i));
    } // for

    return complexResult;
}

std::vector<float> Filter::hammingWindow(int n) {
    std::vector<float> window;

    for (int i = 0; i < n; i++) {
        window.push_back(0.54 - 0.46 * cos(2 * M_PI * i / n));
    }
    
    return window;
}

std::vector<std::complex<float>> Filter::complexHammingWindow(int n) {
    std::vector<std::complex<float>> window;

    for (int i = 0; i < n; i++) {
        std::complex<float> val = (0.54 - 0.46 * cos(2 * M_PI * i / n));
        window.push_back(val);
    }
    
    return window;
}

std::vector<std::complex<float>> Filter::applyFilterByConv(std::vector<std::complex<float>> seq, std::vector<std::complex<float>> fil, int seqLen, int filLen) {
    std::vector<std::complex<float>> output;
    // loop over each val in ouput
    for (int n = 0; n < (seqLen + filLen); n++) {
        std::complex<float> sum = 0;
        // loop over each val in filter
        for (int i = 0; i < filLen; i++) {
            // check filter and sequence array bounds
            if (((n-i) >= 0) && ((n-i) < seqLen)) {
                sum += (fil[i] * seq[n-i]);
            } // else     
        } // for
        output.push_back(sum);
    } // for

    return output;
}

double * Filter::applyFilterByConv(double * seq, double * fil, int seqLen, int filLen) {
    double * output = (double*) malloc((seqLen + filLen) * sizeof(double*));
    // loop over each val in ouput
    for (int n = 0; n < (seqLen + filLen); n++) {
        double sum = 0;
        // loop over each val in filter
        for (int i = 0; i < filLen; i++) {
            // check filter and sequence array bounds
            if (((n-i) >= 0) && ((n-i) < seqLen)) {
                sum += ((*(fil + i)) * (*(seq + n - i)));
            } // else     
        } // for
        *(output + n) = sum;
    } // for

    return output;
}

std::vector<float> Filter::getPrinciple(std::vector<std::vector<float>> inp, std::vector<float> arg) {
    std::vector<float> res;

    for (int i = 0; i < inp.size(); i++) {
        float max = inp[i][0];
        int argmax = 0;
        for (int j = 1; j < inp[i].size(); j++) {
            if (inp[i][j] > max) {
                max = inp[i][j];
                argmax = j;
            }
        } // for
        res.push_back(arg[argmax]);
    } // for

    return res;
}