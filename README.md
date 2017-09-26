# MDLJ
Molecular Dynamics Simulation of particles interacting via the Lennard Jones Potential.

C++ source file to be found in directory "src".

Functions:
Initialize - sets up atoms in fcc crystal lattice with random initial velocities.
Force - calculates total force of configuration.
record - saves a snapshot of configuration to file.
ave_var -calculates the running averages and variance of quantities.
update_gr - accumulates values for the radial distribution function.
order_parameter - calculates the value of the positional order parameter.
vacf - calculates the velocity autocorrelation function.  

To run:
In the directory with lj.cpp, create new folders "1", "0.7".

Run as :
g++ lj.cpp -o lj
./lj


Results and Plots:

Resulting files and plots and scripts for generating the plots to be found in the directory "files"

To plot, run scripts as:

gnuplot <script-filename>

