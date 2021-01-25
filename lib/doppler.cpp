#include <iostream>
#include <cmath>
#include <unordered_map>
#include "doppler.h"

using namespace std;

# define C 343

Doppler::Doppler(int transFreq)
    : transmitFreq(transFreq) 
{ }

Doppler::~Doppler() { }

std::vector<float> Doppler::measureProjectileVelocity(std::vector<float> receivedFreq) {
    std::vector<float> estimProjVel;

    for (int i = 0; i < receivedFreq.size(); i++) {
        float vel = C * (Doppler::transmitFreq - receivedFreq[i])/(Doppler::transmitFreq + receivedFreq[i]);
        estimProjVel.push_back(vel);
    }
    
    return estimProjVel;
}

std::vector<float> Doppler::estimDistanceTravelled(std::vector<float> estimProjVel, std::vector<float> timeSteps) {
    std::vector<float> estimProjDis;
    float dis = 0;
    estimProjDis.push_back(dis);

    for (int i = 1; i < estimProjVel.size(); i++) {
        dis += (timeSteps[i]-timeSteps[i-1])*estimProjVel[i];
        estimProjDis.push_back(dis);
    }
    
    return estimProjDis;
}

std::vector<float> Doppler::measureTransverseDistance(std::vector<float> estimProjVel, std::vector<float> timeSteps) {
    std::vector<float> transDis;

    // here we break causality by estimating the steady
    // state object velocity by finding the mode
    float velSS = estimProjVel[0];
    unordered_map<float, int> modeMap;

    for (int i = 0; i < estimProjVel.size(); i++) {
        modeMap[estimProjVel[i]]++;
    } // for

    int maxCount = 0;
    unordered_map<float, int>::iterator it;

    for (it = modeMap.begin(); it != modeMap.end(); it++) {
        if (it->second > maxCount) {
            maxCount = it->second;
            velSS = it->first;
        } // if
    } // for

    // hypotenuse distance
    std::vector<float> hypDis = Doppler::estimDistanceTravelled(estimProjVel, timeSteps);

    // estimate cos(angle) between reciever and object
    std::vector<float> alpha;

    for (int i = 0; i < estimProjVel.size(); i++) {
        alpha.push_back(acos(estimProjVel[i] / velSS));
    } // for


    for (int i = 0; i < hypDis.size(); i++) {
        // prevent NANs
        float cosAlpha = cos(alpha[i] - M_PI/2);
        if (alpha[i] <= 0 || alpha[i] > M_PI || isnan(alpha[i])) {
            transDis.push_back(0);
        } else {
            transDis.push_back(hypDis[i] * cosAlpha);
        } // else
    } // for
    

    return transDis;
}