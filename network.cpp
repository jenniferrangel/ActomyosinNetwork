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

    cout << "My number of filaments is " << get_num_filaments() << endl;
}

//=======================================================
// Getters and Setters
//=======================================================

//=======================================================
// Destructor
//=======================================================
Network::~Network() {
	cout << "Destructor executed!" << endl;
}
//=========================
//End of network.cpp

