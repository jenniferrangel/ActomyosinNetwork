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
//***Update node positions via Langevin Equation***//
void Node::update_Position(int Ti){
    cout << "New force exerted on node: " << new_total_force << endl;
    cout << "Position: " << my_location << endl;
    cout << "Drag coefficient: " << drag_coeff << endl;
    cout << "dt: " << dt << endl;

    if(isnan(new_total_force.get_X()) || isnan(new_total_force.get_Y())){
        cout << "New force is NaN in update_Position() Ti = " << Ti << endl;
    } else if(isnan(dt)){
        cout << "dt is NaN in update_Position() Ti = " << Ti << endl;
    } else if(isnan(drag_coeff)){
        cout << "drag coeff is NaN in update_Position() Ti = " << Ti << endl;
    }

    my_location += new_total_force*(dt/drag_coeff);

    cout << "Updated node position = " << my_location << endl;

    return;
}

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

//***Calculate force functions***//
//=================================
//***calculates total force on the actin node
void Actin_Node::calculate_Forces(int Ti){
    cout << "Calculating total force..." << endl;

    //------Linear force from left and right neighbors------
    //since 1st node doesn't have L nhbr => F^L = 0. Since last node doesn't have R nhbr => F^R = 0.
    //else F^spring = F^L + F^R
    cout << "Linear spring force being calculated..." << endl;
    Coord F_linear = calc_Linear_Force();

    //------Bending force from left and right neighbors------

    //------Stochastic/random force------


    //total force actin on node
    new_total_force = F_linear;

    return;
}

//***Linear spring force of left and right neighbors
Coord Actin_Node::calc_Linear_Force(){
    Coord F_lin;
    
    cout << "Calc linear from LEFT neighbor" << endl;
    Coord F_left = linear_Spring_Equation(left_neighbor);
    cout << "F_left = " << F_left << endl;

    cout << "Calc linear from RIGHT neighbor" << endl;
    Coord F_right = linear_Spring_Equation(right_neighbor);
    cout << "F_right = " << F_right << endl;

    F_lin = F_left + F_right;
    cout << "F_linear = " << F_lin << endl;

    return F_lin;
}

Coord Actin_Node::linear_Spring_Equation(shared_ptr<Actin_Node> node){
    if(node == NULL){
        cout << "ERROR: Accessing a NULL pointer. Will abort!!!" << endl;
        exit(1);
    }

    Coord F_linear;

    //compute the left/right vector and length
    Coord diff_vec = node->get_Node_Location() - my_location;
    double diff_length = diff_vec.length();
    cout << "diff_vec = " << diff_vec << endl;
    cout << "diff_length = " << diff_length << endl;

    //when we don't have a L or R nbhr, diff_vec = (0,0) and diff_length = 0 so F_linear = (0,0)
    if(diff_length == 0){
        cout << "Avoiding 0 in the denominator. Returning (0,0)" << endl;
        return Coord(0,0);
    }

    //F_spring = k(|node_i - node_j| - l_equi)((node_i - node_j)/|node_i - node_j|)
    F_linear = (diff_vec/diff_length)*k_linear_actin*(diff_length - this->actin_spring_equi_len);
    cout << "Linear Force: " << F_linear << endl;

    return F_linear;

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