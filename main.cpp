#ifndef _MAIN_INCLUDE_
#define _MAIN_INCLUDE_

// Main.cpp
//============================

//===========================
// Include Dependencies
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include <stdio.h>
#include <memory>
#include <random>
#include <functional>
#include <chrono>
#include "params.h"
#include "coord.h"
#include "node.h"
#include "filament.h"
#include "network.h"
//#include "rand.h"
//==============================

using namespace std;

//FREQUENTLY CHANGED VALUES

//Actin parameters:
//============================================
//parameters for linear springs:
double ACTIN_SPRING_EQUI_LEN = 0.2;     //(microns) actin linear spring equilibrium length
double K_LINEAR_STIFF_ACTIN = 150.0;     //(pN) The actin linear spring coefficient
//parameters for bending springs:
double K_BEND_STIFF_ACTIN = 0.0414;     //(pN*micron^2) The actin bending coefficientat every node triplet
double THETA_EQUI_ANGLE = 3.14159265;  // The equilibrium angle which is pi
double ACTIN_DRAG_COEFF = 0.216;        //(pN s/micron)
double VISCOSITY_CYTOPLASM = 0.301;     //(Pa s OR pN s/micron)
double ACTIN_PERSISTENCE_LEN = 10.0;    //(microns)    //(microns)

//Myosin parameters:
//============================================
double F_MYO_PULLING = 4.0;        //(pN) myosin pulling force
double MYO_DRAG_COEFF = 400.0;     //(pN s/micron) myosin drag coefficient

//Actomyosin parameters:
//============================================
double ACTIN_MYO_CONNECT_RADIUS = 0.15;  //(microns) actin and myosin need to be within this distance to form a connection
double K_LINEAR_STIFF_ACTIN_MYO_CONN = 2.0;  //(pN/micron) spring constant for the actin-myosin connection
double ACTIN_MYO_CONN_SPRING_EQUI_LEN = 0.0;  //(microns) equilibrium length of actin-myosin connection (linear spring)

double kB = 1.38064852e-5;              //Boltzmann Const. (micron pN/K)
double TEMPERATURE = 300.0;             //Kelvin
double dt = 0.0005;                     //time step (seconds)



int main(){

    //cout << "Hello I am in main" << endl;

    //start the random number generator
    std::random_device rd;
	std::mt19937 gen(rd());
	srand(time(0));
    
    string init_network = "actomyo_initial_condition.txt";

    Network actomyosin_Network(init_network,gen);


    return 0;

}

#endif