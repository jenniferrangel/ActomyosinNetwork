//node.cpp
//
//=========================
#include <iostream>
#include <vector>
#include <iterator>
#include <cmath>
#include <memory>
#include <algorithm>
//=========================
#include "params.h"
#include "coord.h"
#include "node.h"
#include "filament.h"
#include "network.h"
#include "externs.h"
//=========================

//==========================================================
/** class Node Functions **/
//==========================================================

//==============================
// Constructors
//==============================
Node::Node(Coord location){
    my_location = location;
    new_total_force = Coord();  //sets force to be 0 (0,0)
}

//==============================
//Getters and setters:
//==============================
void Node::set_Drag_Coeff(double drag){
    this->drag_coeff = drag;
    return;
}

void Node::set_My_Node_Rank(int node_num){
    this->my_node_rank = node_num;
    return;
}

//==============================
//Functions
//==============================

//==============================
//Destructor:
//==============================
Node::~Node() {
    cout << "Node destructor executed!" << endl;
}

//===========================================

//==========================================================
/** class Actin Node Functions **/
//==========================================================

//==============================
// Constructors
//==============================
Actin_Node::Actin_Node(Coord location, shared_ptr<Filament> my_filament) : Node(location){
    this->my_filament = my_filament;
}

Actin_Node::Actin_Node(Coord location, shared_ptr<Filament> my_filament, shared_ptr<Actin_Node> left_nbh, shared_ptr<Actin_Node> right_nbh) : Node(location){
    this->my_filament = my_filament;
    this->left_neighbor = left_nbh;
    this->right_neighbor = right_nbh;

    //compute the angle this node is associated with
    calc_Current_Angle();
}

//==============================
//Getters and setters:
//==============================
void Actin_Node::set_My_Filament(shared_ptr<Filament> filament){
    this->my_filament = filament;
    return;
}

void Actin_Node::set_Left_Neighbor(shared_ptr<Actin_Node> left_nbh){
    this->left_neighbor = left_nbh;
    return;
}

void Actin_Node::set_Right_Neighbor(shared_ptr<Actin_Node> right_nbh){
    this->right_neighbor = right_nbh;
    return;
}

void Actin_Node::set_K_Linear_Actin(double k_linear){
    this->k_linear_actin = k_linear;
    return;
}

void Actin_Node::set_Actin_Equi_Len(double equi_len){
    this->actin_spring_equi_len = equi_len;
    return;
}

void Actin_Node::set_K_Bend_Actin(double k_bend){
    this->k_bend_actin = k_bend;
    return;
}

void Actin_Node::calc_Current_Angle(){
    //Recall that cos(theta) = (r_left . r_right)/(|r_left| |r_right|) and theta = acos((r_left . r_right)/(|r_left| |r_right|))
    //NOTE: the FIRST and LAST nodes in a filament don't have an angle i.e. theta = 0
        //* since the first node does not have a left neighbor and last node does not have a right neighbor, we set the neighbors equal to the nodes
        //* first node location = first node's left neighbor location
        //* last node location = last node's right neighbor location
        //* this will give us an angle = 0 for first and last position

    //First compute the left and right vectors
    Coord left_vec = get_Left_Neighbor()->get_Node_Location() - get_Node_Location();
    Coord right_vec = get_Right_Neighbor()->get_Node_Location() - get_Node_Location();

    //Then compute the lengths
    double left_length = left_vec.length();
    double right_length = right_vec.length();

    //avoid a zero in the denominator or an empty
    if(left_length * right_length == 0){
        my_current_angle = 0;
        cout << "Zero in denominator! Overlapping nodes in cal_Current_Angle()! My angle: " << my_current_angle << endl;
        return;
    }else if (isnan(left_length * right_length)){
        cout << "NaN in cal_Current_Angle() Filament " << get_My_Filament()->get_Filament_Num() << endl;
		exit(1);
    }

    double costheta = left_vec.dot(right_vec) / (left_length * right_length);
    double theta = acos( min(max(costheta,-1.0), 1.0) );

    //acos only returns angles between 0 and pi, but the angle can be > pi. So we need to calculate the cross product to get correct angle
        // if crossproduct > 0 then angle is between 0 and pi
        // if crossproduct < 0 then angle is > pi and we need to do 2pi-theta to get correct angle
    
    double crossProduct = left_vec.cross(right_vec);

    if(crossProduct < 0.0){
        theta = 2 * pi - theta;
    }

    //update the angle and cross product in protected member variables
    my_current_angle = theta;
    cout << "My angle: " << theta << endl;
    cross_product = crossProduct;

    return;  
}

void Actin_Node::set_Current_Angle(double my_angle){
    this->my_current_angle = my_angle;
    return;
}

void Actin_Node::set_Equi_Angle(double new_angle){
    this->equi_angle = new_angle;
    return;
}

//==============================
//Functions:
//==============================
 void Actin_Node::sound_Off_Node_Info(){
    cout << "Parent filament #: " << get_My_Filament()->get_Filament_Num() << endl;
    cout << "Node #: " << get_My_Node_Rank() << endl;
    cout << "Location: " << get_Node_Location() << endl;
    cout << "Left Neighbor: " << get_Left_Neighbor()->get_Node_Location() << endl;
    cout << "Right Neighbor: " << get_Right_Neighbor()->get_Node_Location() << endl;
    cout << "Drag coeff: " << get_Drag_Coeff() << endl;
    cout << "Linear spring coeff: " << get_K_Linear_Actin() << endl;
    cout << "Linear equilibrium length: " << get_Actin_Equi_Len() << endl;
    cout << "Bending spring coeff: " << get_K_Bend_Actin() << endl;
    cout << "Current angle: " << get_Current_Angle() << endl;
    cout << "Equilibrium angle: " << get_Equi_Angle() << endl;
    cout << "Total force: " << get_Total_Force() << endl;

    return;
 }

//==============================
//Destructor:
//==============================
Actin_Node::~Actin_Node(){
    cout << "Actin node destructor called!" << endl;
}

//==========================================================
/** class Myosin Node Functions **/
//==========================================================

//==============================
// Constructors
//==============================

//==============================
//Getters and setters:
//==============================

//==============================
//Destructor:
//==============================

//==========================================================
// End of node.cpp