#pragma once 
#include <iostream>
#include <utility> 

class AbstractWalker{
    private: 
        int seed = 5;
    protected: 
        //Updates the particle position one block in a random direction step
        void updatePositionRandomwalk(double setpos[], double pos[], int val);

        //Abstract method that defines the jump length for the generating process. May be unit length, 
        //or  may be random drawn from some kind of distribution
        virtual int generateJumpLength();

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
        std::pair<double, double> updateGaussianJump(double mu, double sigma);

    public: 
        const double mu;
        const double sigma;

        GaussianWalker(double mu = 0.0, double sigma = 1.0) : mu(mu), sigma(sigma) {}

};



int main(void){
    GaussianWalker foo = GaussianWalker(0.0, 4.0);
    std::cout << "hello world" << std::endl;
    return 0; 
}