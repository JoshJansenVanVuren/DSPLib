#ifndef DOPPLER_H
#define DOPPLER_H

#include <vector>


class Doppler {
    public:
        /**
         * @brief Construct a new Doppler object
         * 
         * @param transFreq Known transmission frequency of Doppler radar
         */
        Doppler(int transFreq);
        ~Doppler();

        /**
         * @brief Measures the velocity (m/s) of a projectile from a stationary source
         * 
         * @param receivedFreq Time series of recieved frequencies
         * @return std::vector<float> 
         */
        std::vector<float> measureProjectileVelocity(std::vector<float> receivedFreq);
        
        /**
         * @brief Estimates distance travelled by discretely
         * integrating a time series of object velocities
         * 
         * @param estimProjVel Time series of objects velocity
         * @param timeSteps Time series in seconds (s)
         * @return std::vector<float> 
         */
        std::vector<float> estimDistanceTravelled(std::vector<float> estimProjVel, std::vector<float> timeSteps);

        /**
         * @brief Estimates transverse distance between moving object
         * and a stationary observer by using the cosine error
         * 
         * @param estimProjVel Time series of objects velocity
         * @param timeSteps Time series in seconds (s)
         * @return std::vector<float> 
         */
        std::vector<float> measureTransverseDistance(std::vector<float> estimProjVel, std::vector<float> timeSteps);
    
    private:
        int transmitFreq;
};

#endif // STFT_H