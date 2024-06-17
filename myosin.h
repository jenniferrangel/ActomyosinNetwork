//myosin.h
//===================
// Inlcude Guards
#ifndef _MYOSIN_H_INCLUDED_
#define _MYOSIN_H_INCLUDED_
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
#include "filament.h"
#include <boost/random.hpp>
//===================
// Myosin Mini-filament Class Declaration
class Myosin: public enable_shared_from_this<Myosin> {
    private:
        Network* my_network;

        int myosinNumber; //myosin rank
        double myosin_drag_coeff;

        int num_myosin_nodes;
        vector<shared_ptr<Myosin_Node>> myosin_nodes;

    public:
    //Constructor:
    //=======================
    //Constructor for the case when the myosin nodes are predefined
    Myosin(Network* network, int myosin_Number);

    //make nodes with predefined positions
    void make_myosin_nodes(vector<Coord> nodes);
       
    //Destructor
    //=======================
    ~Myosin();


    //Getters & setters
    //=======================
    Network* get_Network(){return my_network;}

    //get/set the myosin node rank (i.e. is it the first or second node)
    int get_Myosin_Num(){return myosinNumber;}
    void set_Myosin_Num(const int rank);

    //get/set damping/drag coefficient for myosin mini-filaments
    double get_Myosin_Drag_Coeff(){return myosin_drag_coeff;}
    void set_Myosin_Drag_Coeff(double new_drag);

    //get/set number of myosin mini-filament nodes:
    int get_Num_Myosin_Nodes(){return num_myosin_nodes;}
    void set_Num_Myosin_Nodes(int number_myosin_nodes);

    //get myosin nodes
    void get_Myosin_Nodes_Vec(vector<shared_ptr<Myosin_Node>>& myosins);


    //Functions
    //=======================

    //***Calculate forces


    //***Update node positions via Langevin eqn


    //***Functions for VTK output



    //***Functions for data output
    //This function prints out the myosin node locations
    void print_Myosin_Locations(ofstream& ofs, int Ti);
    void print_Myosin_Node_Data(ofstream& ofs, int Ti);
    void print_MiniFilament_Data(ofstream& ofs, int Ti);


};


// End Myosin Mini-filament Class
//===================

#endif
