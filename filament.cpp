//filament.cpp:
//===================
// Forward Declarations
// Include Dependencies
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <random>
#include <memory>
#include "params.h"
#include "coord.h"
#include "node.h"
//#include "rand.h"
#include "filament.h"
#include "network.h"
#include <boost/random.hpp>
//===================

//==========================================================
/** Filament Class Member Functions**/

//==============================
// Constructors
//==============================
Filament::Filament(Network* network, int fil_Number, bool firstIsBarbed, bool lastIsBarbed, Coord initNode){
    this->my_network = network;
    this->filamentNumber = fil_Number;
    this->firstNodeIsBarbed = firstIsBarbed;
    this->lastNodeIsBarbed = lastIsBarbed;
    this->firstNode = initNode;

    //set the damping/drag coefficient
    this->actin_drag_coeff = ACTIN_DRAG_COEFF;
    cout << "My drag coeff is: " << actin_drag_coeff << endl;

    num_actin_nodes = 0;
    set_Initial_Num_Actin_Nodes(INIT_NUM_ACTIN_NODES);
    cout << "My initial # actin nodes is: " << get_Initial_Num_Actin_Nodes() << endl;

}

void Filament::make_nodes(){
    cout << "NEED TO FINISH make_nodes() function!!!" << endl;
}

//==============================
//Getters and setters:
//==============================
void Filament::set_Filament_Num(const int rank){
    this->filamentNumber = rank;
    return;
}

void Filament::set_First_Node_Polarity(bool polarity){
    this->firstNodeIsBarbed = polarity;
    return;
}

void Filament::set_Last_Node_Polarity(bool polarity){
    this->lastNodeIsBarbed = polarity;
    return;
}

void Filament::set_Actin_Drag_Coeff(double new_drag){
    this->actin_drag_coeff = new_drag;
    return;
}

void Filament::set_Initial_Num_Actin_Nodes(int initNumber){
    this->initial_num_actin_nodes = initNumber;
    return;
}

void Filament::set_Num_Actin_Nodes(int number_actin_nodes){
    this->num_actin_nodes = number_actin_nodes;
    return;
}

void Filament::get_Actin_Nodes_Vec(vector<shared_ptr<Actin_Node>>& actins){
    actins = this->actin_nodes;
    return;
}

void Filament::add_Actin_Node_Vec(shared_ptr<Actin_Node> curr){
    this->actin_nodes.push_back(curr);
    this->num_actin_nodes++;
    return;
}

//==============================
//Functions:
//==============================

//==============================
//Destructor:
//==============================
Filament::~Filament(){
    cout << "Filament destructor executed!" << endl;
}


//==========================================================
// End of filament.cpp