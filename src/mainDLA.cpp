#include <GL/glut.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <string>
#include <string.h>
#include <time.h>

#include "../headers/DLASystem.h"
#include "../headers/Window.h" 
#include <iostream>

using namespace std;

// functions which are needed for openGL go into a namespace so that we can identify them
namespace drawFuncs {
  void handleKeypress(unsigned char key, int x, int y);
  void display(void);
  void update(int val);
  void introMessage();
}

// this is a global pointer, which is how we access the system itself
DLASystem *sys;


// this is just a help message
void drawFuncs::introMessage() {
        sys->setRunning();
        glutTimerFunc(0, drawFuncs::update, 0);
        sys->setFast();
        cout << "Keys (while in graphics window):" << endl << "  q or e to quit (or exit)" << endl;
        cout << "  h to print this message (help)" << endl;
        cout << "  u for a single update" << endl;
        cout << "  g to start running (go)" << endl;
        cout << "  p to pause running" << endl;
        cout << "  s to run in slow-mode" << endl;
        cout << "  f to run in fast-mode" << endl;
        cout << "  r to clear everything (reset)" << endl;
        cout << "  z to pause and zoom in" << endl;
        cout << "  w or b to change background colour to white or black" << endl;
        cout << "" << endl;
        cout << "" << endl;
}

// openGL function deals with the keyboard
void drawFuncs::handleKeypress(unsigned char key, int x, int y) {
	switch (key) {
  case 'h':
    drawFuncs::introMessage();
    break;
	case 'q':
	case 'e':
		cout << "Exiting..." << endl;
    // delete the system
    delete sys;
		exit(0);
		break;
  case 'p':
    cout << "pause" << endl;
    sys->pauseRunning();
    break;
  case 'g':
    cout << "go" << endl;
    sys->setRunning();
    glutTimerFunc(0, drawFuncs::update, 0);
    break;
  case 's':
    cout << "slow" << endl;
    sys->setSlow();
    break;
  case 'w':
    cout << "white" << endl;
    sys->setWinBackgroundWhite();
    break;
  case 'b':
    cout << "black" << endl;
    sys->setWinBackgroundBlack();
    break;
  case 'f':
    cout << "fast" << endl;
    sys->setFast();
    break;
  case 'r':
    cout << "reset" << endl;
    sys->Reset();
    break;
  case 'z':
    cout << "zoom" << endl;
    sys->pauseRunning();
    sys->viewAddCircle();
    break;
  case 'u':
    cout << "upd" << endl;
    //sys->Update(logfile);
    break;
  case 'o':
    sys->LogRadius();
    break;
  case 'l':
    int random_run = sys->rgen.randomInt(1000000000);
    ofstream logfile("output" + to_string(random_run) + ".txt");
    for(auto var : sys->LogfileRows){
      logfile << var << endl;
    }
    logfile.close();
    break;
	}
  // tell openGL to redraw the window
	glutPostRedisplay();
}

// this function gets called whenever the algorithm should do its update
void drawFuncs::update(int val) {
  int wait;  // time to wait between updates (milliseconds)
  
  if ( sys->running ) {
    if ( sys->slowNotFast == 1)
      wait = 10;
    else
      wait = 0;
    
    sys->Update();
    
    // tell openGL to call this funtion again after "wait" milliseconds
    glutTimerFunc(wait, drawFuncs::update, 0);
  }
  
}

// this function redraws the window when necessary
void drawFuncs::display() {
    	//  Clear the window or more specifically the frame buffer...
      //  This happens by replacing all the contents of the frame
      //  buffer by the clear color (black in our case)
      glClear(GL_COLOR_BUFFER_BIT);

      // this puts the camera at the origin (not sure why) with (I think) z axis out of page and y axis up
      // there is also the question of the GL perspective which is not set up in any clear way at the moment
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      gluLookAt(0.0, 0.0, 1.0,   /* camera position */
        0.0, 0.0, -1.0,        /* point to look at */
        0.0, 1.0, 0.0);		   /* up direction */

      //sys->DrawSpheres();
      sys->DrawSquares();

      //  Swap contents of backward and forward frame buffers
      glutSwapBuffers();
}

int main(int argc, char **argv) {
  // turn on glut
	glutInit(&argc, argv);
  int window_size[] = { 480,480 };
  string window_title("simple DLA simulation");

  // create a window
  Window *win = new Window(window_size,window_title);
  //Check that seed is actually passed to program
	if (argc == 1){
		throw std::invalid_argument( "No command line args were passed. Please pass an int for an RNG seed for the first argument" );
	};
  // this is the seed for the random numbers
  int seed = std::stoi(argv[1]);
  cout << "setting seed " << seed << endl;

  // create the system
  sys = new DLASystem(win, seed);
  
  // print the "help" message to the console
  drawFuncs::introMessage();

  // tell openGL how to redraw the screen and respond to the keyboard
	glutDisplayFunc( drawFuncs::display );

	glutKeyboardFunc(drawFuncs::handleKeypress );
  
  // tell openGL to do its first update after waiting 10ms
  int wait = 10;
  int val = 0;


	glutTimerFunc(wait, drawFuncs::update, val);
  // start the openGL stuff

 	glutMainLoop();
  delete sys;
  return 0;

}
