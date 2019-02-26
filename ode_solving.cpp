// ode_solving.cpp

// Differential equations describing populations in an apartment building of 500 people infested by zombies:
// dS/dt = - BSZ - ESK
// dK/dt = - CKZ + ESK
// dZ/dt = BSZ + CKZ - AKZ

// Outputs a netcdf file that contains (S, K, Z) population vs. time for two values of the initial Z0 values:
// one in which the zombies win, and one in which they all disappear (Z<1). 


#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <boost/array.hpp>
#include <vector>
#include <netcdf>
#include <rarray>
#include <rarrayio>
using namespace boost::numeric::odeint;
using namespace std;
using namespace netCDF;


// defining constants and initial conditions
const double A = 0.03;  // rate at which zombies are killed by K
const double B = 0.02;  // rate at which regular people are turned into Z
const double C = 0.01;  // rate as which zombie killers are turned into Zombies
const double E = 0.015; // rate at which zombie killers teach regular people how to kill zombies
const double K0 = 9; // number of people who can kill zombies

int i = 0; // integration steps
int saves = 0; // counter for number of saves
double num = 50; // max time
int max_observes = 200; // max number of observes we're doing

// defining 3D array to hold data history
// dimensions: tSKZ, max 200 observes, 2 saves
rarray<double,3> history(4, max_observes, 2);

// defining boost array to hold x = ( S, K, Z )
typedef boost::array<double, 3> state_type;

// saves data in array
void add_to_array(const state_type &x , const double t){
    history[0][i][saves] = t; history[1][i][saves] = x[0]; history[2][i][saves] = x[1]; history[3][i][saves] = x[2]; 
    i = i+1;
}

// x = ( S, K, Z )
void zombie_odes( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = -(B*x[0]*x[2])-(E*x[0]*x[1]);
    dxdt[1] = -(C*x[1]*x[2])+(E*x[0]*x[1]);
    dxdt[2] = -(B*x[0]*x[2])+(C*x[1]*x[2])-(A*x[1]*x[2]);
}

// Function writes netCDF file containing 3D matrix
int netCDF_write(rarray<double,3>& array_to_print) {
   // Finding the length of each dimension
   int nx = array_to_print.extent(0);
   int ny = array_to_print.extent(1);
   int nz = array_to_print.extent(2);
   double dataOut[nx][ny][nz]; // Place to hold data
    
   // Copy data
   for(int x = 0; x < nx; x++){
        for(int y = 0; y < ny; y++){
             for(int z = 0; z < nz; z++){
                 dataOut[x][y][z] = array_to_print[x][y][z];
             }
        }  
   }
    
   // Create the netCDF file.
   NcFile dataFile("output.nc", NcFile::replace);
   // Create the three dimensions.
   NcDim xDim = dataFile.addDim("x",nx);
   NcDim yDim = dataFile.addDim("y",ny);
   NcDim zDim = dataFile.addDim("z",nz);
    
   // Define the variable
   vector<NcDim> dims(3);
   dims[0] = xDim;
   dims[1] = yDim;
   dims[2] = zDim;
   // Create the data variable.
   NcVar data = dataFile.addVar("data", ncDouble, dims);
   // Put the data in the file.
   data.putVar(&dataOut); // writing all the data in one operation
    
   // add an attribute
   dataFile.putAtt("Creation date:", "26 Feb 2019");
   return 0; 
}

// Function reads and formats netCDF file containing 3D matrix according to the zombie ODE solver
int netCDF_read() {
    // Specify the netCDF file. 
    NcFile dataFile("output.nc", NcFile::read);
    
    // Read the three dimensions.
    NcDim xDim = dataFile.getDim("x");
    NcDim yDim = dataFile.getDim("y");
    NcDim zDim = dataFile.getDim("z");
    int nx = xDim.getSize(); 
    int ny = yDim.getSize(); 
    int nz = zDim.getSize(); 

    // Retrieve the variable named "data"
    NcVar data = dataFile.getVar("data");
    // Put the data in a var.
    double dataOut[nx][ny][nz];
    data.getVar(&dataOut);
    
    // Formatting the printout
    for(int z=0; z<nz; z++){
        if(z<1){
            cout << "With 16 initial zombies, they lose" << endl;
        }
        else{
            cout << endl; cout << endl;
            cout << "With 263 initial zombies, they win" << endl;
        }
        cout << "time \t S \t K \t Z" << endl;
        
        for (int y = 0; y < ny; y++) {
            // want to break out of the loop if the data is all 0's
            if ((dataOut[1][y][z]==0)&&(dataOut[2][y][z]==0)&&(dataOut[3][y][z]==0)){
                break;
            }
            // if the data is not all 0's, print it!
            else{
            cout << dataOut[0][y][z] << '\t' << dataOut[1][y][z] << '\t' << dataOut[2][y][z] << '\t' << dataOut[3][y][z] << endl;
            }
        }
    }
    return 0; 
}


int main() { 
  // First case:
  double Z0 = 16; // number of zombies: will be varying this value {16}
  double S0 = 491 - Z0; // number of regular people who can't kill zombies
  state_type x = {S0, K0, Z0}; // initial conditions
  // integrate needs (system, x0, t0, t1, dt, observer)
  integrate( zombie_odes , x , 0.0 , num , 0.01 ,  add_to_array);
    
  // Second case:
  i = 0; // resetting integration steps
  Z0 = 263; // new number of zombies {263}
  S0 = 491 - Z0; // number of regular people who can't kill zombies
  x = {S0, K0, Z0}; // new initial conditions
  saves = 1; // counter for number of saves
  integrate(zombie_odes , x , 0.0 , num , 0.01 , add_to_array);
    
  netCDF_write(history);
  netCDF_read();
}
 
