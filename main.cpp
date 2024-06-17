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
#include "myosin.h"
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
double MYOSIN_SPRING_EQUI_LEN = 0.15; //NEED TO FIND ACTUAL VALUE
double K_LINEAR_STIFF_MYOSIN = 600.0; //NEED TO FIND VALUE SO DOESNT STRETCH

//Actomyosin parameters:
//============================================
double ACTIN_MYO_CONNECT_RADIUS = 0.15;  //(microns) actin and myosin need to be within this distance to form a connection
double K_LINEAR_STIFF_ACTIN_MYO_CONN = 2.0;  //(pN/micron) spring constant for the actin-myosin connection
double ACTIN_MYO_CONN_SPRING_EQUI_LEN = 0.0;  //(microns) equilibrium length of actin-myosin connection (linear spring)

double kB = 1.38064852e-5;              //Boltzmann Const. (micron pN/K)
double TEMPERATURE = 300.0;             //Kelvin
double dt = 0.0005;                     //time step (seconds)

bool PREDEFINED_INITIAL_NODES = true;
int PRINT_VTKS = true; // ./batchGenerator -par -PRINT <1 or 0>
int NUM_STEPS_PER_FRAME = 1;


int main(int argc, char* argv[]){

    //Read in the names of folders where we store outputs
    //*****************************************************
    //Reads in name of folder that stores vtk files. This is "Animation" folder.
    string animation_folder = argv[1];

    //Reads in the name of folder that stores data output. 
    //This is DataOutput/Locations
    string locations_folder = argv[2];
    //Folders to hold the Node, Filament and Network level data
    string node_data_folder = argv[3];
    string filament_data_folder = argv[4];
    string myosin_data_folder = argv[5];
    string network_data_folder = argv[6];

    //start the random number generator
    //*****************************************************
    std::random_device rd;    //gives a random seed from operating system
	std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
	srand(time(0));          //ensures that the seed value is different each time the program is run
    
    //instantiate the network: making the filaments and their respective nodes
    //*****************************************************
    cout << "Generating network" << endl;

    //***if we want to generate nodes within code
    // string init_network = "actomyo_initial_condition.txt";
    // Network actomyosin_Network(init_network,gen,PREDEFINED_INITIAL_NODES);

    //***if the initial condition file already has the position of all nodes:
    //input file that already has all of the nodes predefined and PREDEFINED_INITIAL_NODES=true
    //string init_network_with_nodes = "actomyo_nodes_initial.txt";
    string init_network_with_nodes = "actomyosin_net_nodes_initial.txt";
    Network actomyosin_Network(init_network_with_nodes,gen,PREDEFINED_INITIAL_NODES);

    cout << "Finished creating the network" << endl;

    //for double-checking purposes
    actomyosin_Network.sound_Off_All_Node_Info();
    actomyosin_Network.sound_Off_Neighbors();

    //Variables for outputs
    //*****************************************************
    //variables for vtk files
    //========================
    int digits;
    string Filename;
    string initial = "/Actin_Filament";
    string Number;
    string format = ".vtk";	
	ofstream ofs_anim;
	int out = 0;

    //variables for dataoutput
    //========================
    //For actin filaments:
    //int digits2;
	string Locations; //this is the filename
    string locations_initial = "/Actin_Locations_";
    //string Number2;
	ofstream ofs_loc;
	int out2 = 0;

    string Node_Data;
    string node_data_initial = "/Actin_Node_Data_";
    ofstream ofs_node_data;
    int out3 = 0;

    string Filament_Data;
    string filament_data_initial = "/Filament_Data_";
    ofstream ofs_filament_data;
    int out4 = 0;

    //For myosin mini-filaments:
    string Locations_Myosin; //this is the filename
    string locations_initial_myo = "/Myosin_Locations_";
    //string Number2;
	ofstream ofs_loc_myosin;
	int out6 = 0;

    string Node_Data_Myosin;
    string node_data_myosin_initial = "/Myosin_Node_Data_";
    ofstream ofs_node_data_myosin;
    int out7 = 0;

    string Myosin_Data;
    string myosin_data_initial = "/Myosin_MiniFilament_Data_";
    ofstream ofs_myosin_data;
    int out8 = 0;

    //For actomyosin network:
    string Network_Data;
    string network_data_initial = "/Network_Data_";
    ofstream ofs_network_data;
    int out5 = 0;

    //Start the loop
    //*****************************************************
    int Ti = 0;
    int final_time = 4;

    while(Ti < final_time){
        cout << "Entered while loop" << endl;

        //if node positions are updated before vtk and data output files printed, we won't get the info for the initial condition

        //Calculate the forces on filament nodes
        //*****************************************************
        /* cout << "Calculating forces..." << endl;
        actomyosin_Network.calculate_New_Forces(Ti);
        cout << "Forces calculated!" << endl; */

        // update filament node positions
        //*****************************************************
/*         cout << "Updating locations..." << endl;
        actomyosin_Network.update_Positions(Ti);
        cout << "Locations updated!!" << endl;

        //For double checking purposes
        actomyosin_Network.sound_Off_Neighbors(); */
        

        //Print vtk files and to DataOutput
        //*****************************************************
        //Printing out the VTK files
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

                cout << "Opening the file in the animation folder" << endl;
                ofs_anim.open(Filename.c_str());
                actomyosin_Network.print_VTK_File(ofs_anim);
                ofs_anim.close();

                out++;
            }
        }
        cout << "VTK printed" << endl;

        //Data output from simulations
        if(Ti % NUM_STEPS_PER_FRAME == 0){
            //Printing actin filament data info:
            //===================================
            //Printing out the node locations to DataOutput
            Locations = locations_folder + locations_initial + to_string(out2) + ".txt";
            ofs_loc.open(Locations.c_str());
            cout << "Locations output file opened..." << endl;
            actomyosin_Network.locations_Output(ofs_loc, Ti);
            ofs_loc.close();
            out2++;

            //Printing out the node data 
            Node_Data = node_data_folder + node_data_initial + to_string(out3) + ".txt";
            ofs_node_data.open(Node_Data.c_str());
            cout << "Node data output file opened..." << endl;
            actomyosin_Network.node_Data_Output(ofs_node_data, Ti);
            ofs_node_data.close();
            out3++;

            //Printing out the filament data
            Filament_Data = filament_data_folder + filament_data_initial + to_string(out4) + ".txt";
            ofs_filament_data.open(Filament_Data.c_str());
            cout << "Filament data output file opened..." << endl;
            actomyosin_Network.filament_Data_Output(ofs_filament_data, Ti);
            ofs_filament_data.close();
            out4++;

            //Printing myosin mini-filament data info:
            //===================================
            //Printing out the myosin node locations to DataOutput
            Locations_Myosin = locations_folder + locations_initial_myo + to_string(out6) + ".txt";
            ofs_loc_myosin.open(Locations_Myosin.c_str());
            cout << "Myosin Locations output file opened..." << endl;
            actomyosin_Network.locations_Myosin_Output(ofs_loc_myosin, Ti);
            ofs_loc_myosin.close();
            out6++;

            //Printing out the myosin node data 
            Node_Data_Myosin = node_data_folder + node_data_myosin_initial + to_string(out7) + ".txt";
            ofs_node_data_myosin.open(Node_Data_Myosin.c_str());
            cout << "Myosin Node data output file opened..." << endl;
            actomyosin_Network.myosin_Node_Data_Output(ofs_node_data_myosin, Ti);
            ofs_node_data_myosin.close();
            out7++;

            //Printing out the myosin mini-filament data
            Myosin_Data = myosin_data_folder + myosin_data_initial + to_string(out8) + ".txt";
            ofs_myosin_data.open(Myosin_Data.c_str());
            cout << "Myosin Mini-Filament data output file opened..." << endl;
            actomyosin_Network.Myosin_Minifilament_Data_Output(ofs_myosin_data, Ti);
            ofs_myosin_data.close();
            out8++;

            //Printing out the network data
            //===============================
            Network_Data = network_data_folder + network_data_initial + to_string(out5) + ".txt";
            ofs_network_data.open(Network_Data.c_str());
            cout << "Network data output file opened..." << endl;
            actomyosin_Network.network_Data_Output(ofs_network_data, Ti);
            ofs_network_data.close();
            out5++;

        }

        //This allows file '0' i.e. at time = 0 to have the info for the initial condition

        //Calculate the forces on filament nodes
        //*****************************************************
        cout << "Calculating forces..." << endl;
        actomyosin_Network.calculate_New_Forces(Ti);
        cout << "Forces calculated!" << endl;

        // update filament node positions
        //*****************************************************
        cout << "Updating locations..." << endl;
        actomyosin_Network.update_Positions(Ti);
        cout << "Locations updated!!" << endl;

        //For double checking purposes
        actomyosin_Network.sound_Off_Neighbors();
        

        Ti++;

    }
	
	
    return 0;

}

#endif