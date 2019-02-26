// ode_solving.cpp

// Differential equations describing populations in an apartment building of 500 people infested by zombies:
// dS/dt = - BSZ - ESK
// dK/dt = - CKZ + ESK
// dZ/dt = BSZ + CKZ - AKZ

// Want to numerically solve the set of ODEs for several values of Z0 to determine the number of Z0 above which
// the zombies will win. 
// Outputs a netcdf file that contains (S, K, Z) population vs. time for two values of the initial Z0 values:
// one in which the zombies win (everyone is zombies), and one in which they all disappear (Z<1). 


#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <boost/array.hpp>
#include <rarray>
#include <rarrayio>
using namespace boost::numeric::odeint;
using namespace std;


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


int main() { 
  // First case:
  double Z0 = 16; // number of zombies: will be varying this value
  double S0 = 491 - Z0; // number of regular people who can't kill zombies
  state_type x = {S0, K0, Z0}; // initial conditions
  // integrate needs (system, x0, t0, t1, dt, observer)
  integrate( zombie_odes , x , 0.0 , num , 0.01 ,  add_to_array);
    
  // Second case:
  i = 0; // resetting integration steps
  Z0 = 263; // new number of zombies 
  S0 = 491 - Z0; // number of regular people who can't kill zombies
  x = {S0, K0, Z0}; // new initial conditions
  saves = 1; // counter for number of saves
  integrate(zombie_odes , x , 0.0 , num , 0.01 , add_to_array);
}
