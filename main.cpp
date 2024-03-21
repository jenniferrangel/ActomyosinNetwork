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

int PRINT_VTKS = true; // ./batchGenerator -par -PRINT <1 or 0>
int NUM_STEPS_PER_FRAME = 1;



int main(int argc, char* argv[]){

    //Read in the names of folders where we store outputs
    //*****************************************************
    //Reads in name of folder that stores vtk files. This is "Animation" folder.
    string animation_folder = argv[1];

    //Reads in the name of folder that stores data output. This is DataOutput/Locations
    string locations_folder = argv[2];

    //start the random number generator
    //*****************************************************
    std::random_device rd;
	std::mt19937 gen(rd());
	srand(time(0));
    
    string init_network = "actomyo_initial_condition.txt";

    //instantiate the network: making the filaments and their respective nodes
    //*****************************************************
    cout << "Generating network" << endl;
    Network actomyosin_Network(init_network,gen);
    cout << "Finished creating the network" << endl;

    //for double-checking purposes
    actomyosin_Network.sound_Off_All_Node_Info();

    //Variables for outputs
    //*****************************************************
    //variables for vtk files
    int digits;
    string Filename;
    string initial = "/Actin_Filament";
    string Number;
    string format = ".vtk";	
	ofstream ofs_anim;
	int out = 0;

    //variables for dataoutput

    //Start the loop
    //*****************************************************
    int Ti = 0;
    int final_time = 3;

    while(Ti < final_time){
        cout << "Entered while loop" << endl;

        if(PRINT_VTKS){

            if(Ti % NUM_STEPS_PER_FRAME == 0){
                digits = ceil(log10(out + 1));
				if (digits == 1 || digits == 0) {
					Number = "0000" + to_string(out);
				}	
				else if (digits == 2) {
					Number = "000" + to_string(out);
				}	
				else if (digits == 3) {
					Number = "00" + to_string(out);
				}
				else if (digits == 4) {
					Number = "0" + to_string(out);
				}

                Filename = animation_folder + initial + Number + format;

                cout << "Opening the animation folder" << endl;
                ofs_anim.open(Filename.c_str());
                actomyosin_Network.print_VTK_File(ofs_anim);
                ofs_anim.close();

                out++;
            }
        }

        Ti++;

    }
	
	
    return 0;

}

#endif