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

//for the case when nodes are predefined:
Filament::Filament(Network* network, int fil_Number, bool firstIsBarbed, bool lastIsBarbed){
    this->my_network = network;
    this->filamentNumber = fil_Number;
    this->firstNodeIsBarbed = firstIsBarbed;
    this->lastNodeIsBarbed = lastIsBarbed;
    
    //set the damping/drag coefficient
    this->actin_drag_coeff = ACTIN_DRAG_COEFF;
    cout << "My drag coeff is: " << actin_drag_coeff << endl;

    num_actin_nodes = 0;
}

void Filament::make_nodes(){
    cout << "Making the actin nodes!" << endl;

    int num_Initial_Actin_Nodes = INIT_NUM_ACTIN_NODES;
    double increment = ACTIN_SPRING_EQUI_LEN;

    //Make all of the actin nodes
    double current_x;
    double current_y;
    Coord location = this->firstNode;

    current_x = firstNode.get_X();
    current_y = firstNode.get_Y();
    location = Coord(current_x, current_y);

    //make the first node
    shared_ptr<Filament> this_filament = shared_from_this();
    shared_ptr<Actin_Node> previousActin = make_shared<Actin_Node>(location, this_filament);
    shared_ptr<Actin_Node> currentActin = previousActin;
    actin_nodes.push_back(previousActin); //add it to vector of nodes
    num_actin_nodes++;    //increment number of nodes

    shared_ptr<Actin_Node> original(previousActin); //to save the starter node
    //to fix the neighbors of first and last node
    shared_ptr<Actin_Node> first(previousActin);
    shared_ptr<Actin_Node> last = previousActin;

    //make the successive actin nodes
    for(int i = 0; i < num_Initial_Actin_Nodes - 1; i++){
      current_x = firstNode.get_X() + increment*(i+1);
      current_y = firstNode.get_Y();
      location = Coord(current_x, current_y);

      shared_ptr<Actin_Node> new_node = make_shared<Actin_Node>(location, this_filament);
      currentActin = new_node;
      actin_nodes.push_back(currentActin);
      num_actin_nodes++;

      //setting the left and right neighbors
      previousActin->set_Right_Neighbor(currentActin);
      currentActin->set_Left_Neighbor(previousActin);

      previousActin = currentActin;
      last = currentActin;
    }

    //The first and last nodes do not have a pair of neighbors
    //The FIRST node does NOT have a LEFT neighbor so we will set the left neighbor to itself. Ex if node is (0,0) then left neighbor = (0,0)
    original->set_Left_Neighbor(first);
    
    //The LAST node does NOT have a RIGHT neighbor so we will set the right neighbor to itself. Ex if node is (1,0) then right neighbor = (1,0)
    currentActin->set_Right_Neighbor(last);

    //Now we set the private node member variables
    vector<shared_ptr<Actin_Node>> actins; //will store actin nodes here
    this->get_Actin_Nodes_Vec(actins);
    //the drag coefficient is inherited from the filament
    double new_drag_coeff = this->get_Actin_Drag_Coeff();
    // double spring_equi_len = 0.0;
    // double k_linear = 0.0;
    // double k_bend = 0.0;
    int rank = 0;

    for(unsigned int i = 0; i < actins.size(); i++){
        rank = i;
        actins.at(i)->set_Drag_Coeff(new_drag_coeff);
        actins.at(i)->set_Actin_Equi_Len(ACTIN_SPRING_EQUI_LEN);
        actins.at(i)->set_K_Linear_Actin(K_LINEAR_STIFF_ACTIN);
        actins.at(i)->set_K_Bend_Actin(K_BEND_STIFF_ACTIN);
        actins.at(i)->set_Equi_Angle(THETA_EQUI_ANGLE);

        actins.at(i)->set_My_Node_Rank(rank);
    }

    //update the actin filament angles
    update_Actin_Angles();

    //set the barbed end
    set_Barbed_End();

    cout << "Finished making the actin nodes!" << endl;

    return;
}

//make nodes with predefined positions
void Filament::make_nodes(vector<Coord> nodes){
    cout << "Making the actin nodes!" << endl;

    double current_x;
    double current_y;
    Coord location;
    Coord node;
    shared_ptr<Filament> this_filament = shared_from_this();
    shared_ptr<Actin_Node> previousActin;
    shared_ptr<Actin_Node> currentActin;

    shared_ptr<Actin_Node> original; //to save the starter node
    //to fix the neighbors of first and last node
    shared_ptr<Actin_Node> first;
    shared_ptr<Actin_Node> last;

    //make the successive actin nodes
    for(unsigned int i = 0; i < nodes.size(); i++){
        if(i==0){
            this->firstNode = nodes.at(i);
            current_x = firstNode.get_X();
            current_y = firstNode.get_Y();
            location = Coord(current_x, current_y);
            shared_ptr<Actin_Node> new_node = make_shared<Actin_Node>(location, this_filament);
            original = new_node;
            first = new_node;
            previousActin = new_node;
            actin_nodes.push_back(previousActin);
            num_actin_nodes++;
        }else{
            node = nodes.at(i);
            current_x = node.get_X();
            current_y = node.get_Y();
            location = Coord(current_x, current_y);
            shared_ptr<Actin_Node> new_node = make_shared<Actin_Node>(location, this_filament);
            currentActin = new_node;
            actin_nodes.push_back(currentActin);
            num_actin_nodes++;

            //setting the left and right neighbors
            previousActin->set_Right_Neighbor(currentActin);
            currentActin->set_Left_Neighbor(previousActin);

            previousActin = currentActin;
            last = currentActin;
        }

    }

    //The first and last nodes do not have a pair of neighbors
    //The FIRST node does NOT have a LEFT neighbor so we will set the left neighbor to itself. Ex if node is (0,0) then left neighbor = (0,0)
    original->set_Left_Neighbor(first);
    
    //The LAST node does NOT have a RIGHT neighbor so we will set the right neighbor to itself. Ex if node is (1,0) then right neighbor = (1,0)
    currentActin->set_Right_Neighbor(last);

    //Now we set the private node member variables
    vector<shared_ptr<Actin_Node>> actins; //will store actin nodes here
    this->get_Actin_Nodes_Vec(actins);
    //the drag coefficient is inherited from the filament
    double new_drag_coeff = this->get_Actin_Drag_Coeff();
    // double spring_equi_len = 0.0;
    // double k_linear = 0.0;
    // double k_bend = 0.0;
    int rank = 0;

    for(unsigned int i = 0; i < actins.size(); i++){
        rank = i;
        actins.at(i)->set_Drag_Coeff(new_drag_coeff);
        actins.at(i)->set_Actin_Equi_Len(ACTIN_SPRING_EQUI_LEN);
        actins.at(i)->set_K_Linear_Actin(K_LINEAR_STIFF_ACTIN);
        actins.at(i)->set_K_Bend_Actin(K_BEND_STIFF_ACTIN);
        actins.at(i)->set_Equi_Angle(THETA_EQUI_ANGLE);

        actins.at(i)->set_My_Node_Rank(rank);
    }

    //update the actin filament angles
    update_Actin_Angles();

    //set the barbed end
    set_Barbed_End();

    cout << "Finished making the actin nodes!" << endl;

    return;

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

void Filament::set_Barbed_End(){
    //get polarity of first and last node
    bool firstNodeBarbed = this->get_First_Node_Polarity();
    bool lastNodeBarbed = this->get_Last_Node_Polarity();

    //access the nodes in the filament
    vector<shared_ptr<Actin_Node>> actins; //will store actin nodes here
    this->get_Actin_Nodes_Vec(actins);

    Coord firstNode;
    Coord lastNode;

    if(firstNodeBarbed){
        cout << "First node is the barbed end" << endl;

        firstNode = actins.at(0)->get_Node_Location();
        this->barbedEnd = firstNode;
        cout << "Barbed end = " << get_Barbed_End() << endl;
    }
    else if(lastNodeBarbed){
        cout << "Last node is the barbed end" << endl;

        lastNode = actins.at(actins.size() - 1)->get_Node_Location();
        this->barbedEnd = lastNode;
        cout << "Barbed end = " << get_Barbed_End() << endl;
    }

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
void Filament::update_Actin_Angles(){
    vector<shared_ptr<Actin_Node>> actins; //will store actin nodes here
    this->get_Actin_Nodes_Vec(actins);

//#pragma omp parallel for schedule(static,1)
    for(unsigned int i = 0; i < actins.size(); i++){
        actins.at(i)->calc_Current_Angle();
    }
    cout << "Finished calculating actin filament angles" << endl;
    return;
}

//***Calculate forces***//
void Filament::calculate_New_Forces(int Ti){
    vector<shared_ptr<Actin_Node>> actins; //will store actin nodes here
    this->get_Actin_Nodes_Vec(actins);

    cout << "Calculating forces..." << endl;

    //#pragma omp parallel
       // {
    //#pragma omp for schedule(static,1)
            for(unsigned int i = 0; i < actins.size(); i++){
                actins.at(i)->calculate_Forces(Ti);
            }
        //}

    return;
}

//***Update node positions via Langevin Equation***//
void Filament::update_Node_Positions(int Ti){
    vector<shared_ptr<Actin_Node>> actins; //will store actin nodes here
    this->get_Actin_Nodes_Vec(actins);

    //#pragma omp parallel 
	    //{	
    //#pragma omp for schedule(static,1)
            for(unsigned int i = 0; i < actins.size(); i++){
                cout << "Updating locations..." << endl;
                actins.at(i)->update_Position(Ti);
            }
        //}

    update_Actin_Angles();
    
    return;

}

//***Functions for VTK output****//
void Filament::print_VTK_Points(ofstream& ofs, int& count){
    for(unsigned int i = 0; i < actin_nodes.size(); i++){
        Coord location = actin_nodes.at(i)->get_Node_Location();
        ofs << location.get_X() << ' ' << location.get_Y() << ' ' << 0 << endl;
        count++;
    }

    return;
}

void Filament::print_VTK_BarbedEnds(ofstream&ofs){
    Coord barbed = get_Barbed_End();
    ofs << barbed.get_X() << ' ' << barbed.get_Y() << ' ' << 0 << endl;

    return;

}

//***Functions for data output***//
//This function prints out the node locations
void Filament::print_Location(ofstream& ofs, int Ti){
    //ofs <<"Filament#" << ' ' << "Node#" << ' ' << "Location (x, y, z)" << endl;
    for(unsigned int i = 0; i < actin_nodes.size(); i++){
        Coord location = actin_nodes.at(i)->get_Node_Location();
        ofs << this->get_Filament_Num() << ' ' << actin_nodes.at(i)->get_My_Node_Rank() << ' '<< location.get_X() << ' ' << location.get_Y() << ' ' << 0 << endl;

    }
    return;
}

//This function prints out the node data
void Filament::print_Node_Data(ofstream& ofs, int Ti){
    //ofs <<"Filament#" << ' ' << "Node#" << ' ' << "Location (x, y, z)" << ' ' <<  "Left Nbr (x, y, z)" << ' ' <<  "Right Nbr (x, y, z)" << ' ' << "Drag Coeff" << ' ' << "k_linear" << ' ' << "Equi Len" << ' ' << "k_bend" << ' ' << "theta" << ' ' << "theta_equi" << endl;
    for(unsigned int i = 0; i < actin_nodes.size(); i++){
        Coord location = actin_nodes.at(i)->get_Node_Location();
        Coord left_neighbor = actin_nodes.at(i)->get_Left_Neighbor()->get_Node_Location();
        Coord right_neighbor = actin_nodes.at(i)->get_Right_Neighbor()->get_Node_Location();

        ofs << this->get_Filament_Num() << ' ' << actin_nodes.at(i)->get_My_Node_Rank() << ' '<< location.get_X() << ' ' << location.get_Y() << ' ' << 0 << ' ' << left_neighbor.get_X() << ' ' << left_neighbor.get_Y() << ' ' << 0 << ' ' << right_neighbor.get_X() << ' ' << right_neighbor.get_Y() << ' ' << 0 << ' ' << actin_nodes.at(i)->get_Drag_Coeff() << actin_nodes.at(i)->get_K_Linear_Actin() << ' ' << actin_nodes.at(i)->get_Actin_Equi_Len() << ' ' << actin_nodes.at(i)->get_K_Bend_Actin() << ' ' << actin_nodes.at(i)->get_Current_Angle() << ' ' << actin_nodes.at(i)->get_Equi_Angle() << endl;
    }

    return;
}

//This function prints out the filament data
void Filament::print_Filament_Data(ofstream& ofs, int Ti){
    ofs << this->get_Filament_Num() << ' ' << this->get_Num_Actin_Nodes() << ' ' << this->get_First_Node_Polarity() << ' ' << this->get_Last_Node_Polarity() << endl;

    return;
}

//==============================
//Destructor:
//==============================
Filament::~Filament(){
    cout << "Filament destructor executed!" << endl;
}


//==========================================================
// End of filament.cpp