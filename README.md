# 3D-GP-MPIOMP
3D-GP-MPIOMP is a FORTRAN project designed to numerically solve the Gross-Pitaevskii equation (GPE) in three dimensions (3D) using MPI and OpenMP.

Solving the GPE allows for qualitatively accurate simulations of Bose-Einstein Condensates (BECs) at zero temperature.  
3D-GP solves the GPE using 4th order Runge-Kutta time stepping on a grid of points with regular spacing.
User defined grid spacing and time step is supported.

### Requirements
3D-GP requires NetCDF and NetCDF-Fortran installations to run. NetCDF is used to save compressed data files.
On Ubuntu you can install the required packages by running:
```
sudo apt-get install libnetcdf11 libnetcdff6 libnetcdf-dev libnetcdff-dev
```

You can also create a local installation of NetCDF by downloading the source files from the
[NetCDF website](https://www.unidata.ucar.edu/software/netcdf/) and compiling them yourself.
If you do this, be sure to build with parallel NetCDF support and make a note of the installation
location as you will need it to install 3D-GP.

### Installation
* Clone the git repository somewhere.

* Run `./install` to setup. You can also run  `./install <install-dir>` to install to any desired installation location.

* New terminals should now be able to run  `make3dgp` anywhere.

### Running a Simulation
* Create a new simulation directory
* Enter the directory and run `make3dgp` to set up a simulation at your current location.
* Edit `run.sh` to set the number of parallel processes.
* Type `./run.sh` to start the simulation.
* The simulation status is printed to screen.

### Editing Parameters
To run a simulation with custom parameters
*  Create a simulation as in **Running a Simulation** but do not run `./run.sh`.
*  Edit `params.in` to include any parameters you wish to change from their defaults and run `make3dgp` again to reflect the changes.
*  Run `./run.sh` to start the simulation.
