//
//  DLASystem.cpp
//
#include <math.h>
#include <string.h>
#include "../headers/DLASystem.h"
#include<cmath>
#include<algorithm>
#include <cmath>
#include <chrono>
#include <thread>

#define protonMass 1.6726219e-27
// colors
namespace colours {
	GLfloat blue[] = { 0.1, 0.3, 0.9, 1.0 };   // blue
	GLfloat red[] = { 1.0, 0.2, 0.1, 0.2 };   // red
	GLfloat green[] = { 0.3, 0.6, 0.3, 1.0 };     // green
	GLfloat paleGrey[] = { 0.7, 0.7, 0.7, 1.0 };     // green
	GLfloat darkGrey[] = { 0.2, 0.2, 0.2, 1.0 };     // green
}

// this function gets called every step,
//   if there is an active particle then it gets moved,
//   if not then add a particle
void DLASystem::Update() {

	if (lastParticleIsActive == 1) {
		moveLastParticle();
	}

	//At end of simulation write to output text file
	else if(numParticles == endNum){		
		string fp = "data/force_vector/force_vector_0.5/output" + to_string(seed) + ".txt";
		std::cout << "writing simulation data for seed " << to_string(seed) << " to" << fp << std::endl;
		ofstream logfile(fp);
		for(auto var : LogfileRows){
		logfile << var << endl;
		}
		logfile.close();
		exit(0);
	}

	else if (numParticles < endNum) {
		addParticleOnAddCircle();
		setParticleActive();
	}
	if (lastParticleIsActive == 0 || slowNotFast == 1)
		glutPostRedisplay(); //Tell GLUT that the display has changed

}


std::vector<std::pair<int, int>> DLASystem::generateRing(std::pair<int, int> particlePosition, int radius, int** grid) {
    std::vector<std::pair<int, int>> coordinatePositiveRing; // initialize number of points in the ring to zero
	int x = particlePosition.first;
	int y = particlePosition.second; 

    // loop over a square region of width 2r+1 centered at (x, y)
    for (int i = x - radius; i <= x + radius; i++) {
        for (int j = y - radius; j <= y + radius; j++) {
            // check if (i, j) is within the ring
            if (std::sqrt((i - x)*(i - x) + (j - y)*(j - y)) < radius) {
                // check if (i, j) is within the bounds of the grid
                if (i >= 0 && i < 1600 && j >= 0 && j < 1600) {
                    //Check if the array position is 1
					if(grid[i][j] == 1){
						coordinatePositiveRing.emplace_back(i , j);
					}
                }
            }
        }
	}
	return coordinatePositiveRing;
}


std::pair<int,int> DLASystem::isPathClear(int x1, int y1, int x2, int y2) {
	//cout << "path coords " << x1 << "," << y1 << "," << x2 << "," << y2 << endl;
    // Calculate the x and y distances between the points
    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);

    // Determine the direction of travel in each dimension
    int dirX = (x2 > x1) ? 1 : -1;
    int dirY = (y2 > y1) ? 1 : -1;

    // Calculate the error in each dimension
    int error = dx - dy;

    // Track the current position along the path
    int x = x1;
    int y = y1;

    // Traverse along the path from (x1,y1) to (x2,y2)
    while (x != x2 || y != y2) {
        // Check if the current position contains a particle
        //if (x+dirX == 0 && y+dirY == 0){
		//if(grid[x+800][y+800] == 1){
		//	return make_pair(x,y);
        //}
        // Update the error and move in the x direction
        int error2 = error * 2;
        if (error2 > -dy) {
            error -= dy;
            x += dirX;
        }
        // Update the error and move in the y direction
        if (error2 < dx) {
            error += dx;
            y += dirY;
        }
    }
    // No particles found along the path
    return std::make_pair(x,y);
}


//Given the magnitude of the force, find how many jumps this means the random walker should move
//This acts as a discretisation of the force magnitude 
int findJumpSize(int forceVectorMagnitude){
	vector<pair<int, int>> forceCutOffs = {{24/protonMass,6}, {20/protonMass,5}, {16/protonMass,4}, {12/protonMass,3}, {8/protonMass,2}, {4/protonMass,1}};
	for (auto& forceVal : forceCutOffs){
		if (signbit(forceVal.first - abs(forceVectorMagnitude)) ==1){
			if (signbit(forceVectorMagnitude) == 1){return forceVal.second * -1;}
			return forceVal.second;
		}
	}
	if (signbit(forceVectorMagnitude) == 1){return-1;}
	else{return 1;}

}


//Given a particle position and the mean force vector that is acting on that particle from an attractive 
//potential, calculate the jump size that needs to be walked in each direction to find the new position, 
//Find the geodesic path and try to move there. If a particle (in the cluster) is in the way, collide with it and 
//cease walking along the geodesic line 
std::pair<int, int> DLASystem::findNewPosition(pair<int, int> particlePosition, pair<double, double> meanForceVector){
	std::pair<int, int> jumpVector = {findJumpSize(meanForceVector.first), findJumpSize(meanForceVector.second)};
	std::pair<int, int> DesirednewParticlePosition = {particlePosition.first + jumpVector.first, particlePosition.second + jumpVector.second};
	
  	std::pair<int, int> GeodesicParticlePositon = isPathClear(particlePosition.first, particlePosition.second, DesirednewParticlePosition.first, DesirednewParticlePosition.second);
	if (GeodesicParticlePositon != DesirednewParticlePosition){
		cout << "collision!" << endl;
	}
	return GeodesicParticlePositon;
}

//Find the yukawa potential at any given position. Define this as a ratio of the magnitude of the brownian (Langevin force)
//This does not check that a particle is in the local area for there to be a force, that has to be done seperately. 
std::pair<double, double> DLASystem::findForceVector(std::pair<int,int> randomwalkerPosition, std::pair<int,int> fractalParticlePosition){
	//if xDiff is positive, force is to right, if Ydiff is positive then force is also up. 
	float xDiff =  randomwalkerPosition.first - fractalParticlePosition.first;
	float yDiff =  randomwalkerPosition.second -fractalParticlePosition.second;
	float totalDiff = sqrt(pow(xDiff,2) + pow(yDiff,2));

	//Find the mag of the force as a Yakuma potential 	
	float forceVector = maximumForceScale * exp(-totalDiff * protonMass)/ pow(totalDiff,2);
	//a = F/m
	forceVector  = forceVector / protonMass; 
	float xUnit = xDiff / totalDiff; 
	float yUnit = yDiff / totalDiff; 

	float forcexComponent = xUnit * forceVector ;
	float forceyComponent = yUnit * forceVector ;

	std::pair<float, float> forceCoords = std::make_pair(forcexComponent, forceyComponent);
	return forceCoords;
}

//Given several vectors in the 2D plane acting on the ranom walking particle, find their mean.
std::pair<double, double> DLASystem::findVectorMean(std::vector<std::pair<double,double>> VectorList){
	vector<double> xVals, yVals; 
	for (const auto& vector : VectorList){
		//Particle force is compared to itself in DLA::System generateRing, this has distance 0 and thus inf force. 
		//This should be 0, so do not add to forces to mean
		if (isnan(vector.first) == 0 || isnan(vector.second) == 0){
			xVals.push_back(vector.first);
			yVals.push_back(vector.second);
		}
	}

	if (xVals.size() != 0){
		float xAverage = accumulate(xVals.begin(), xVals.end(), 0.0/ xVals.size());
		float yAverage = accumulate(yVals.begin(), yVals.end(), 0.0/ yVals.size());
		return std::make_pair(xAverage, yAverage);
	}
	else {return std::make_pair(0, 0);}
}

//Reformat the pair of vector components to a vector of size 4 which
//is the RNG weights for random walking
std::vector<double> DLASystem::convertVectorRngProbability(std::pair<double, double> vectorMean){
	double vectorMaxValInitial = maximumForceScale / 10;
	std::vector<double> vectorProb = {1,1,1,1,};
	if(signbit(vectorMean.first) == 0)
		{vectorProb[0] = vectorProb[0]+vectorMean.first;}
	else
		{vectorProb[1] = vectorProb[1] - vectorMean.first;}

	if(signbit(vectorMean.second) == 0)
		{vectorProb[2] = vectorProb[2]+vectorMean.second;}
	else
		{vectorProb[3] = vectorProb[3] - vectorMean.second;}

	return vectorProb;
}


void DLASystem::clearParticles() {
	// delete particles and the particle list
	for (int i = 0; i < numParticles; i++) {
		delete particleList[i];
	}
	particleList.clear();
	numParticles = 0;
}

// remove any existing particles and setup initial condition
void DLASystem::Reset() {
	// stop running
	running = 0;

	clearParticles();

	lastParticleIsActive = 0;

	// set the grid to zero
	for (int i = 0; i < gridSize; i++) {
		for (int j = 0; j < gridSize; j++) {
			grid[i][j] = 0;
		}
	}
	// setup initial condition and parameters
	addCircle = 10;
	killCircle = 2.0*addCircle;
	clusterRadius = 0.0;
	// add a single particle at the origin
	double pos[] = { 0.0, 0.0 };
	addParticle(pos);

	// set the view
	int InitialViewSize = 40;
	setViewSize(InitialViewSize);

}

// set the value of a grid cell for a particular position
// note the position has the initial particle at (0,0)
// but this corresponds to the middle of the grid array ie grid[ halfGrid ][ halfGrid ]
void DLASystem::setGrid(double pos[], int val) {
	int halfGrid = gridSize / 2;
	grid[(int)(pos[0] + halfGrid)][(int)(pos[1] + halfGrid)] = val;
}

// read the grid cell for a given position
int DLASystem::readGrid(double pos[]) {
	int halfGrid = gridSize / 2;
	return grid[(int)(pos[0] + halfGrid)][(int)(pos[1] + halfGrid)];
}

// check if the cluster is big enough and we should stop:
// to be safe, we need the killCircle to be at least 2 less than the edge of the grid
int DLASystem::checkStop() {
	if (killCircle + 2 >= gridSize / 2) {
		pauseRunning();
		cout << "STOP" << endl;
		glutPostRedisplay(); // update display
		return 1;
	}
	else return 0;
}

// add a particle to the system at a specific position
void DLASystem::addParticle(double pos[]) {
	// create a new particle
	Particle * p = new Particle(pos);
	// push_back means "add this to the end of the list"
	particleList.push_back(p);
	numParticles++;

	// pos coordinates should be -gridSize/2 < x < gridSize/2
	setGrid(pos, 1);
}

// add a particle to the system at a random position on the addCircle
// if we hit an occupied site then we do nothing except print a message
// (this should never happen)
void DLASystem::addParticleOnAddCircle() {
	double pos[2];
	double theta = rgen.random01() * 2 * M_PI;
	pos[0] = ceil(addCircle * cos(theta));
	pos[1] = ceil(addCircle * sin(theta));
	if (readGrid(pos) == 0)
		addParticle(pos);
	else
		cout << "FAIL " << pos[0] << " " << pos[1] << endl;
}

// send back the position of a neighbour of a given grid cell
// NOTE: there is no check that the neighbour is inside the grid,
// this has to be done separately...
void DLASystem::setPosNeighbour(double setpos[], double pos[], int val) {
		switch (val) {
		case 0:
			setpos[0] = pos[0] + 1.0;
			setpos[1] = pos[1];
			break;
		case 1:
			setpos[0] = pos[0] - 1.0;
			setpos[1] = pos[1];
			break;
		case 2:
			setpos[0] = pos[0];
			setpos[1] = pos[1] + 1.0;
			break;
		case 3:
			setpos[0] = pos[0];
			setpos[1] = pos[1] - 1.0;
			break;
		}
}


// if the view is smaller than the kill circle then increase the view area (zoom out)
void DLASystem::updateViewSize() {
	double mult = 1.2;
	if (viewSize < 2.0*killCircle) {
		setViewSize(viewSize * mult);
	}
}

// set the view to be the size of the add circle (ie zoom in on the cluster)
void DLASystem::viewAddCircle() {
	setViewSize(2.0*addCircle);  // factor of 2 is to go from radius to diameter
}

// when we add a particle to the cluster, we should update the cluster radius
// and the sizes of the addCircle and the killCircle
void DLASystem::updateClusterRadius(double pos[]) {

	double rr = distanceFromOrigin(pos);
	if (rr > clusterRadius) {
		clusterRadius = rr;
		// this is how big addCircle is supposed to be:
		//   either 20% more than cluster radius, or at least 5 bigger.
		double check = clusterRadius * addRatio;
		if (check < clusterRadius + 5)
			check = clusterRadius + 5;
		// if it is smaller then update everything...
		if (addCircle < check) {
			addCircle = check;
			killCircle = killRatio * addCircle;
			updateViewSize();
		}
		checkStop();
	}
}


// make a random move of the last particle in the particleList
void DLASystem::moveLastParticle() {
	double newpos[2];
	Particle *lastP = particleList[numParticles - 1];
	//std::this_thread::sleep_for(std::chrono::milliseconds(500));

	//Run vanilla random walk mechanics
	if(condition=="vanilla"){
		int rr = rgen.randomInt(4);
		setPosNeighbour(newpos, lastP->pos, rr);
	}
	//Run single step size yukawa potential mechanics 
	if(condition == "force_vector"){
		std::pair<int, int> position = std::make_pair(lastP->pos[0] + 800, lastP->pos[1] + 800);
		std::vector<std::pair<int, int>> nearbyParticles = generateRing(position, 10, grid);
		std::vector<std::pair<double, double>> forceVectors;
		for(int i = 0; i < nearbyParticles.size(); ++i){
			std::pair<double, double> Vector = findForceVector(nearbyParticles[i], position);
			forceVectors.push_back(Vector);
		}

		std::pair<double, double> vectorMean = findVectorMean(forceVectors);		
		std::vector<double> vectorMeanFourD = convertVectorRngProbability(vectorMean);
		int rr = rgen.weightedRandInt(vectorMeanFourD);
		setPosNeighbour(newpos, lastP->pos, rr);
	}

	//run jump process yukawa process
	if(condition == "force_vector_jump"){
		std::pair<int, int> position = std::make_pair(lastP->pos[0] + 800, lastP->pos[1] + 800);
		std::vector<std::pair<int, int>> nearbyParticles = generateRing(position, 10, grid);
		if (nearbyParticles.size() == 0){
				
			int rr = rgen.randomInt(4);
			setPosNeighbour(newpos, lastP->pos, rr);
		}
		else{
			std::vector<std::pair<double, double>> forceVectors;
			for(int i = 0; i < nearbyParticles.size(); ++i){
				std::pair<double, double> Vector = findForceVector(nearbyParticles[i], position);
				forceVectors.push_back(Vector);
			}

			std::pair<double, double> vectorMean = findVectorMean(forceVectors);		
			std::vector<double> vectorMeanFourD = convertVectorRngProbability(vectorMean);
			position.first = position.first - 800; 
			position.second = position.second - 800; 

			std::vector<std::pair<int, int>> shortRangeNearbyParticles = generateRing(position, 3, grid);
			if (shortRangeNearbyParticles.size() > 0){
				std::pair<int, int> newPosition = findNewPosition(position, vectorMean);
				newpos[0] = newPosition.first;
				newpos[1] = newPosition.second;

			}
			else{
				std::vector<double> vectorMeanFourD = convertVectorRngProbability(vectorMean);
				int rr = rgen.weightedRandInt(vectorMeanFourD);
				setPosNeighbour(newpos, lastP->pos, rr);
			}
		}
	}
	if (distanceFromOrigin(newpos) > killCircle) {
		//cout << "#deleting particle" << endl;
		setGrid(lastP->pos, 0);
		particleList.pop_back();  // remove particle from particleList
		numParticles--;
		setParticleInactive();
	}
	// check if destination is empty
	else if (readGrid(newpos) == 0) {
		setGrid(lastP->pos, 0);  // set the old grid site to empty
		// update the position
		particleList[numParticles - 1]->pos[0] = newpos[0];
		particleList[numParticles - 1]->pos[1] = newpos[1];
		setGrid(lastP->pos, 1);  // set the new grid site to be occupied

		// check if we stick
		if (checkStick()){
			//cout << "stick" << endl;
			setParticleInactive();  // make the particle inactive (stuck)
			updateClusterRadius(lastP->pos);  // update the cluster radius, addCircle, etc.
			if (numParticles % 100 == 0 && logfile.is_open()) {
				logfile << numParticles << " " << clusterRadius << endl;
			}
			bool verb = false;
			string LogRow = LogRadius(verb);
			LogfileRows.push_back(LogRow);
		
		}
	}
	else {
		// if we get to here then we are trying to move to an occupied site
		// (this should never happen as long as the sticking probability is 1.0)
		//cout << "reject " << rr << endl;
		//cout << lastP->pos[0] << " " << lastP->pos[1] << endl;
		//cout << newpos[0] << " " << newpos[1] << " " << (int)newpos[0] << endl;
		//printOccupied();
	}
}

// check if the last particle should stick (to a neighbour)
int DLASystem::checkStick(double StickProb) {
	//Stickprob should not be able to be greater than one, it is a probability
	if (StickProb > 1.0){
		throw std::invalid_argument( "StickProb was set greater than one" );
	};
	Particle *lastP = particleList[numParticles - 1];
	int result = 0;
	// loop over neighbours
	if (condition == "vanilla" || condition == "force_vector"){
		for (int i = 0; i < 4; i++) {
			double checkpos[2];
			setPosNeighbour(checkpos, lastP->pos, i);
			// if the neighbour is occupied...
			if (readGrid(checkpos) == 1){
				if (rgen.random01() < StickProb){
					result = 1;
				}
			};
		}
		return result;
	}

}

// constructor
DLASystem::DLASystem(Window *set_win, int seed_, string condition_) {
	cout << "creating system, gridSize " << gridSize << endl;
	win = set_win;
	numParticles = 0;
	endNum = 10000;
	
	//set rng seed 
	seed = seed_;
	setSeed(seed_);
	condition = condition_;

	// allocate memory for the grid, remember to free the memory in destructor
	grid = new int*[gridSize];
	for (int i = 0; i < gridSize; i++) {
		grid[i] = new int[gridSize];
	}
	slowNotFast = 1;
	// reset initial parameters
	Reset();\
	addRatio = 1.3;   // how much bigger the addCircle should be, compared to cluster radius
	killRatio = 1.5;   // how much bigger is the killCircle, compared to the addCircle

	//vector<int> foobar = {1};

	// this opens a logfile, if we want to...
	//logfile.open("opfile.txt");
}

// destructor
DLASystem::~DLASystem() {
	// strictly we should not print inside the destructor but never mind...
	cout << "deleting system" << endl;
	// delete the particles
	clearParticles();
	// delete the grid
	for (int i = 0; i < gridSize; i++)
		delete[] grid[i];
	delete[] grid;

	if (logfile.is_open())
		logfile.close();

}

double DLASystem::FindFractalDimention(int NumParticles, double ClusterRadius){
	//Function to find the Fractal Dimention
	double a = 1; 
	double LogRadialFraction = log(ClusterRadius / a);
	double LogNumParticles = log(NumParticles);
	//fractal_dim = log(N) / log(R)
	return LogNumParticles / LogRadialFraction; 
}

string DLASystem::LogRadius(bool Verbose){
	//Function for logging
	if (Verbose == true){
		cout << "number_particles:" << particleList.size() << ",";
		cout << "cluster_radius:" << clusterRadius << ","; 
		cout << "fractal_Dimention" << FindFractalDimention(particleList.size(), clusterRadius) << endl;
	}
	
	string output = 
		"number_particles:" + to_string(particleList.size()) + "," + 
		"cluster_radius:" + to_string(clusterRadius);	
	return output;
}



// this draws the system
void DLASystem::DrawSquares() {

	// draw the particles
	double halfSize = 0.5;
	for (int p = 0; p < numParticles; p++) {
		double *vec = particleList[p]->pos;
		glPushMatrix();
		if (p == numParticles - 1 && lastParticleIsActive == 1)
			glColor4fv(colours::red);
		else if (p == 0)
			glColor4fv(colours::green);
		else
			glColor4fv(colours::blue);
		glRectd(drawScale*(vec[0] - halfSize),
			drawScale*(vec[1] - halfSize),
			drawScale*(vec[0] + halfSize),
			drawScale*(vec[1] + halfSize));
		glPopMatrix();
	}

	// print some information (at top left)
	// this ostringstream is a way to create a string with numbers and words (similar to cout << ... )
	ostringstream str;
	str << "num " << numParticles << " size " << clusterRadius;
	// print the string
	win->displayString(str, -0.9, 0.9, colours::red);

	// if we are paused then print this (at bottom left)
	if (running == 0) {
		ostringstream pauseStr;
		pauseStr << "paused";
		win->displayString(pauseStr, -0.9, -0.9, colours::red);
	}

}
