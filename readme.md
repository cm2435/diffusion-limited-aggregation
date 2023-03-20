# DLA C++ Simulation

Some undergrad physics work on implimenting some random process simulations, namely the Ising model and DLA in C++. Credits go to the relevent sources. 

## Development 

This project is developed using the following software requirements
    **g++ (Ubuntu 11.3.0-1ubuntu1~22.04) 11.3.0** 
    **OpenGL version string: 4.6 (Compatibility Profile) Mesa 22.2.5** 
    **Cpp 11.3.0** 



## Running this simulation

1. Download Makefile from this section

2. Put Makefile in same directory as DLASystem.cpp (! Remove the .txt extension if present !)

3. Relative to MS Visual Studio code, these files include the following changes: DLASystem.h Window.h mainDLA.cpp

4. Open the terminal in the same folder as these files

5. Type "make" to build project


6. Run the program by typing "./run"
6. a) Run takes in two command line arguments, the seed and the running type
6. b) The two running types for this project are vanilla and force_vector_jump
6. c) Therefore an example command run would be ./run 5 vanilla 
6. d) To multi process this, run bash.sh. This defaults to run with 24 processes in parallel, so it may crash smaller computers

