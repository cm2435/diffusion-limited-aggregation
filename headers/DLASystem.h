#pragma once

#include <GL/glut.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>
#include <string>
#include <sstream>

#include "Window.h"
#include "Particle.h"
#include "rnd.h"
#include <unordered_map>

using namespace std;


class DLASystem {
  private:
  // these are private variables and functions that the user will not see
  
    Window *win;  // window in which the system is running

    //seed for rng
    int seed;

    //maximum force intensity for the yukakwa potential if used.
    int maximumForceScale = 1.5; 
    
    //Check the geodesic path between two points to jump for collisions 
    std::pair<int,int> isPathClear(int x1, int y1, int x2, int y2);

    //Find the particles new position for a jumping process given a force
    std::pair<int, int> findNewPosition(pair<int, int> particlePosition, pair<double, double> meanForceVector);


    //random walk condition type
    string condition;

    // list of particles
    vector<Particle*> particleList;
    int numParticles;

    // delete particles and clear the particle list
    void clearParticles();

    // size of cluster
    double clusterRadius;
    // these are related to the DLA algorithm
    double addCircle;
    double killCircle;
  
    // size of grid
    static const int gridSize = 1600;
    int **grid;  // this will be a 2d array that stores whether each site is occupied
  
    // the window draws only part of the grid, viewSize controls how much...
    double viewSize;
    double drawScale;
  
    
    // output file (not used at the moment)
    ofstream logfile;
  
    // number of particles at which the simulation will stop
    // (the value is set in constructor)
    int endNum;
  
    // the values of these variables are set in the constructor
    double addRatio;    // how much bigger the addCircle should be, compared to cluster radius
    double killRatio;   // how much bigger is the killCircle, compared to the addCircle

    //check the local ring around the random walker to find any particles within radius r 
    std::vector<std::pair<int, int>> generateRing(std::pair<int, int> particlePosition, int radius, int** grid);

    //Calculate the force weighting between a particle in the fractal and the random walker 
    std::pair<double, double> findForceVector(std::pair<int,int> randomwalkerPosition, std::pair<int,int> fractalParticlePosition);

    //Given a series of force vectors find their mean  
    std::pair<double, double> findVectorMean(std::vector<std::pair<double,double>> VectorList);

    //Convert a force 'vector' into a length 4 weight vector for the RNG generation
    std::vector<double> convertVectorRngProbability(std::pair<double, double> vectorMean);
  
    //Velocity vector probability to be incremented at each timestep
	  //std::vector<double> wewwrgh = {1.0, 1.0, 1.0, 1.0};

  public:
  // these are public variables and functions
    // random number generator, class name is rnd, instance is rgen
    rnd rgen;

    //Vector of strings for the log file 
    vector<string> LogfileRows; 

    // update the system: if there is an active particle then move it,
    // else create a new particle (on the adding circle)
    void Update();

    // draw particles as squares
    void DrawSquares();
  
    // is the simulation running (1) or paused (0) ?
    int running;
  
    // slowNotFast is +1 for slow running, 0 for fast
    int slowNotFast;

    // lastParticleIsActive is +1 if there is an active particle in the system, otherwise 0
    int lastParticleIsActive;

    // constructor
    DLASystem(Window *set_win, int seed, string condition);
    // destructor
    ~DLASystem();
  
    // delete all particles and reset
    void Reset();

    // this sets the seed for the random numbers
    void setSeed(int s) { rgen.setSeed(s); }

    // check whether we should stop (eg the cluster has reached the edge of the grid)
    int checkStop();
  
    // stop/start the algorithm
    void setRunning() { if ( checkStop()==0 ) running = 1; }
    void pauseRunning() { running = 0; }
    // set whether it runs fast or slow
    void setSlow() { slowNotFast = 1; }
    void setFast() { slowNotFast = 0; }
    void setSuperFast() { slowNotFast = -1; }

    // set which part of the grid is visible on the screen
    // basically the screen shows co-ordinates -vv < x < vv
    // where vv is the input value
    void setViewSize(double vv) { viewSize = vv; drawScale = 2.0/viewSize; }
  
    // if the killcircle is almost as big as the view then increase the view
    void updateViewSize();

    // set the view to be the approx size of the addCircle
    void viewAddCircle();

    // if pos is outside the cluster radius then set clusterRadius to be the distance to pos.
    void updateClusterRadius( double pos[] );

    // set and read grid entries associated with a given position
    void setGrid(double pos[], int val);
    int readGrid(double pos[]);

    // return the distance of a given point from the origin
    double distanceFromOrigin(double pos[]) {
      return sqrt( pos[0]*pos[0] + pos[1]*pos[1] );
    }

    // set whether there is an active particle in the system or not
    void setParticleActive()   { lastParticleIsActive = 1; }
    void setParticleInactive() { lastParticleIsActive = 0; }

    // add a particle at pos
    void addParticle(double pos[]);
    // add a particle at a random point on the addCircle
    void addParticleOnAddCircle();

    // assign setpos to the position of a neighbour of pos
    // which neighbour we look at is determined by val (=0,1,2,3)
    void setPosNeighbour(double setpos[], double pos[], int val);

    // this attempts to move the last particle in the List to a random neighbour
    // if the neighbour is occupied then nothing happens
    // the function also checks if the moving particle should stick.
    void moveLastParticle();

    // check whether the last particle should stick
    // currently it sticks whenever it touches another particle
    int checkStick(double StickProb = 1);
    //log the current number of patciles and the DLA radius 
    string LogRadius(bool Verbose = true);

    //find the fractal dimention from the number of particles and the radius 
    double FindFractalDimention(int NumParticles, double ClusterRadius); 

    // set the background colour for the window
    // it would be better for an OOP philosophy to make these member functions for the Window class
    // but we are being a bit lazy here
    void setWinBackgroundWhite() { glClearColor(1.0, 1.0, 1.0, 1.0); }
    void setWinBackgroundBlack() { glClearColor(0.0, 0.0, 0.0, 0.0); }
};
