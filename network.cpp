//network.cpp
//=========================

//=========================
//Include Dependencies
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <memory>
#include <algorithm>
#include <random>
#include "params.h"
#include "coord.h"
#include "filament.h"
#include "network.h"
//=========================
// Public Member Functions for network.cpp

//=======================================================
// Constructors
//=======================================================
//network constructor makes new network from .txt file that is read in

Network::Network(string filename, mt19937 gen){
    num_filaments = 0;

    ifstream ifs(filename.c_str());
    
    if(!ifs.is_open()) {
		cout << filename << " is not available! Please check input file name." << endl;
		return;
	}

    stringstream ss;
    string line;
    string temp;
    char garbage;

    //variables to initialize a filament
    int fil_num;
    bool firstIsBarbed;
    bool lastIsBarbed;
    Coord initialNode;
    double x;
    double y;
    Network* my_network = this;

    //parse through the  input file
    while(getline(ifs,line)){
        ss.str(line);

        getline(ss,temp,':');

        if (temp == "FilamentNumber"){
            ss >> fil_num;
        }
        else if(temp == "FirstIsBardedEnd"){
            ss >> firstIsBarbed;            
        }
        else if(temp == "LastIsBarbedEnd"){
            ss >> lastIsBarbed;
        }
        else if(temp == "Node1"){
            ss >> x >> garbage >> y;
            Coord location(x,y);
            initialNode = location;            
        }
        else if(temp == "End_Filament"){
            //input data is used to create new cell. New cell is pushed onto vector that holds all filaments in the network.
            cout << "Making an actin filament" << endl;

            shared_ptr<Filament> curr = make_shared<Filament>(my_network, fil_num, firstIsBarbed, lastIsBarbed, initialNode);
            num_filaments++;
        }
        ss.clear();
    }
    ifs.close();

    cout << "My number of filaments is " << get_Num_Filaments() << endl;
    cout << "My filament # is: " << fil_num << endl;
    if(firstIsBarbed == true){
        cout << "First node is barbed end: " << firstIsBarbed << endl;
    }
    cout << "Last node is barbed end: " << lastIsBarbed << endl;
    cout << "My location is: " << initialNode << endl;
}

//=======================================================
// Getters and Setters
//=======================================================
//returning or updating the number of actin filaments in the network:
void Network::get_Filaments(vector<shared_ptr<Filament>>& filaments){
    filaments = this->filaments;
    return;
}

void Network::update_Num_Filaments(shared_ptr<Filament>& new_Filament){
   num_filaments++;
   filaments.push_back(new_Filament); 
   return;
}

//=======================================================
// Destructor
//=======================================================
Network::~Network() {
	cout << "Network destructor executed!" << endl;
}
//=========================
//End of network.cpp

