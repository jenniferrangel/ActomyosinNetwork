//filament.h
//===================
// Inlcude Guards
#ifndef _FILAMENT_H_INCLUDED_
#define _FILAMENT_H_INCLUDED_
//===================
// forward declarations
class Network;
//===================
// include dependencies
#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <memory>
#include "params.h"
#include "coord.h"
#include "node.h"
//#include "network.h"
#include "externs.h"
#include <boost/random.hpp>
//===================
// Filament Class Declaration

class Filament: public enable_shared_from_this<Filament> {
    private:
        Network* my_network;

        int filamentNumber; //filament rank
        //polarity of filament: if the first node is barbed end, then last node is pointed end and vice versa
        bool firstNodeIsBarbed;
        bool lastNodeIsBarbed;
        Coord barbedEnd;
        double actin_drag_coeff;
        //seed to start the filament
        Coord firstNode;

        int initial_num_actin_nodes;
        int num_actin_nodes;
        vector<shared_ptr<Actin_Node>> actin_nodes;

    public:
    //Constructor:
    //=======================
    //Constructor for when we want to define nodes within the code using position of initial node
    Filament(Network* network, int fil_Number, bool firstIsBarbed, bool lastIsBarbed, Coord initNode);

    //Constructor for the case when nodes are predefined
    Filament(Network* network, int fil_Number, bool firstIsBarbed, bool lastIsBarbed);

    void make_nodes();
    //make nodes with predefined positions
    void make_nodes(vector<Coord> nodes);

    
    //Destructor
    //=======================
    ~Filament();

    //Getters & setters
    //=======================
    Network* get_Network(){return my_network;}
    
    int get_Filament_Num(){return filamentNumber;}
    void set_Filament_Num(const int rank);

    //get/set polarity for the first and last node i.e. are they the barbed or pointed end;
    bool get_First_Node_Polarity(){return firstNodeIsBarbed;}
    void set_First_Node_Polarity(bool polarity);
    bool get_Last_Node_Polarity(){return lastNodeIsBarbed;}
    void set_Last_Node_Polarity(bool polarity);

    //get/set the barbed end
    Coord get_Barbed_End(){return barbedEnd;}
    void set_Barbed_End();

    //get/set damping/drag coefficient for actin filaments
    double get_Actin_Drag_Coeff(){return actin_drag_coeff;}
    void set_Actin_Drag_Coeff(double new_drag);

    //getset the filament seed
    Coord get_First_Node(){return firstNode;}

    //get number of initial actin filament nodes:
    int get_Initial_Num_Actin_Nodes(){return initial_num_actin_nodes;}
    void set_Initial_Num_Actin_Nodes(int initNumber);

    //get number of actin filament nodes:
    int get_Num_Actin_Nodes(){return num_actin_nodes;}
    void set_Num_Actin_Nodes(int number_actin_nodes);

    //get actin nodes
    void get_Actin_Nodes_Vec(vector<shared_ptr<Actin_Node>>& actins);
    
    //NOT USED--- get actin nodes-added so myosins can have access to these nodes
    // const vector<shared_ptr<Actin_Node>>& get_Actin_Nodes(){return actin_nodes;}

    //add a new actin node to the node vector
    void add_Actin_Node_Vec(shared_ptr<Actin_Node> curr);

    //Functions
    //=======================
    void update_Actin_Angles();
    void update_Barbed_End_Pos();

    //***Calculate forces
    void calculate_New_Forces(int Ti);

    //***Update node positions via Langevin eqn
    void update_Node_Positions(int Ti);

    //***Functions for VTK output
    void print_VTK_Points(ofstream& ofs, int& count);
    void update_Node_VTK_Indices(int& id);
    void print_VTK_BarbedEnds(ofstream& ofs);
    void print_VTK_connections(ofstream& ofs);


    //***Functions for data output
    //This function prints out the node locations
    void print_Location(ofstream& ofs, int Ti);
    void print_Node_Data(ofstream& ofs, int Ti);
    void print_Filament_Data(ofstream& ofs, int Ti);

};

// End Filament Class
//===================

#endif
