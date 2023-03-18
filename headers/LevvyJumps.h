#pragma once 
#include <iostream>
#include <utility> 
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>

//THIS CODE IS FOR LEVVY FLIGHTS AND IS NOT USED AS PARK OF THE YUKAWA POTENTIAL 
class AbstractWalker{
    private: 
        std::random_device rd;

        //Abstract method that defines the jump length for the generating process. May be unit length, 
        //or  may be random drawn from some kind of distribution
        virtual int generateJumpLength();

    protected:
        //rng seed
        int seed = 5;
 
        //Updates the particle position one block in a random direction step
        void updatePositionRandomwalk(double setpos[], double pos[], int val);


        //Function that given a 'jump' radius, can choose a random jump angle and find the closest coordinate to jump to 
        std::pair<int, int> setUnitspherePosition(std::pair<int,int> currentPosition, double radius); 

    public: 
        //Abstract method to define how the final simulation step will update the particle position.
        //Can be implimented to not be one block jumps (gaussian or levy walk) and not of random direction
        //(eg : random field)
        virtual void updatePosition(double setpos[], double pos[], int val);
};


class GaussianWalker : AbstractWalker{
    private: 
        std::mt19937 gen;
        std::normal_distribution<double> dist; // mean=mu, standard deviation=sigma
        
        //Generate the jump length from a gaussian distribution
        int generateJumpLength(){double sample = dist(gen); return sample;}

        //Given the jump length and current position of the particle, find a random direction to jump in
        //by drawing circle of radius jumpLength, and update the particle position to be at a random point on the circle 
        std::pair<double, double> updateGaussianJump(std::pair<int, int> currentPosition, int jumpLength);

    protected:

    public: 
        //Mean and variance for a gaussian distribution of params 
        const double mu;
        const double sigma;

    GaussianWalker(double mu = 0.0, double sigma = 1.0) : mu(mu), sigma(sigma), gen(std::random_device{}()), dist(mu, sigma) {}
};

