//node.cpp
//
//=========================
#include <iostream>
#include <vector>
#include <iterator>
#include <cmath>
#include <memory>
#include <algorithm>
#include <random>
//=========================
#include "params.h"
#include "coord.h"
#include "node.h"
#include "filament.h"
#include "network.h"
#include "externs.h"
#include "myosin.h"
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

void Node::update_VTK_Index(int id){
    vtk_index = id;
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

void Actin_Node::sound_Off_Neighbors(){
    cout << "Node #: " << get_My_Node_Rank() << endl;
    cout << "Location: " << get_Node_Location() << endl;
    cout << "Left Neighbor: " << get_Left_Neighbor()->get_Node_Location() << endl;
    cout << "Right Neighbor: " << get_Right_Neighbor()->get_Node_Location() << endl;

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
    cout << "F_linear = " << F_linear << endl;

    //------Bending force from left and right neighbors------
    cout << "Bending force being calculated..." << endl;
    Coord F_bending = calc_Bending_Force();
    cout << "F_bending = " << F_bending << endl;

    //------Stochastic/random force------
    cout <<"Stochastic/random force being calculated..." << endl;
    Coord F_stoch = calc_Stochastic_Force();
    cout << "F_stoch = " << F_stoch << endl;


    //total force actin on node
    new_total_force = F_linear + F_bending + F_stoch;

    cout << "Total force = " << F_linear << " + " << F_bending << " + " << F_stoch << endl;

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

//***Bending force at every node triplet: First and last node do NOT have an angle
Coord Actin_Node::calc_Bending_Force(){
    Coord F_bend_total;
    
    cout << "Calculating bending F^C..." << endl;
    Coord F_bend_center = bending_Force_Equation_Center();   //this is F^C
    cout << "F_center = " << F_bend_center << endl;

    cout << "Calculating bending F^L..." << endl;
    Coord F_bend_left = bending_Force_Equation_Left();       //this is F^L
    cout << "F_left = " << F_bend_left << endl;
   
    cout << "Calculating bending F^R..." << endl;
    Coord F_bend_right = bending_Force_Equation_Right();     //this is F^R
    cout << "F_right = " << F_bend_right << endl;

    cout << "summing the bending forces..." << endl;
    F_bend_total = F_bend_center + F_bend_left + F_bend_right;

    if(cross_product < 0.0){
        F_bend_total = F_bend_total*(-1);
    }

    return F_bend_total;

}

Coord Actin_Node::bending_Force_Equation_Center(){
    Coord F_center;
    double constant_self;

    double epsilon = 0.0001;

    Coord left_nbr = left_neighbor -> get_Node_Location();
    Coord right_nbr = right_neighbor -> get_Node_Location();

    if (left_nbr == my_location) {
        cout << "First node! Has NO angle and left_nbr = my_loc so F^C = 0!!" << endl;
        return F_center;
    }
    else if (right_nbr == my_location) {
        cout << "Last node! Has NO angle and right_nbr = my_loc so F^C = 0!!" << endl;
        return F_center;
    }
    else if (my_current_angle == 0) {
        //Has no angle so does not experience F^C
        cout << "Not a center node because my_angle = 0 so F^C = 0!" << endl;
        return F_center;
    }
    else if (abs(my_current_angle - pi) < epsilon) {
        //takes care of both numerator and denominator being = 0
        cout << "|my_angle - pi| < epsilon! F^C = 0!" << endl;
        return F_center;
    }

    constant_self = k_bend_actin*(my_current_angle - equi_angle)/(sqrt(1-pow(cos(my_current_angle),2)));
    cout << "constant_self = " << constant_self << endl;

    Coord left_vec = left_neighbor -> get_Node_Location() - my_location;
    Coord right_vec = right_neighbor -> get_Node_Location() - my_location;
    double left_length = left_vec.length();
    double right_length = right_vec.length();

    Coord term1 = (left_vec*(-1))/(left_length*right_length);
    Coord term2 = left_vec*cos(my_current_angle)/pow(left_length,2);
    Coord term3 = (right_vec*(-1))/(right_length*left_length);
    Coord term4 = right_vec*cos(my_current_angle)/pow(right_length,2);

    F_center = (term1 + term2 + term3 + term4) * constant_self;
    cout << "Bending center:" << endl;
    cout << "F_center = (" << term1 << " + " << term2 << " + " << term3 << " + " << term4 << ") * " << constant_self << " = " << F_center << endl;

    return F_center;
}

Coord Actin_Node::bending_Force_Equation_Left(){
    Coord F_left;
    double constant_left;
    double epsilon = 0.0001;

    //First deal with the fact that the first node has NO angle and NO left nbr (set to be itself)
    Coord left_nbr = left_neighbor -> get_Node_Location();
    if (left_nbr == my_location){
        cout << "First node! Have NO angle and NO left neighbor so F^L = 0!!" << endl;
        return F_left;
    }

    //The first node has no angle, no left neighbor (set to be itself) and no LL neighbor
    //This is for the 2nd node, where F^L = 0 because it's L nbr (1st node) has NO angle
    Coord left_left_nbr = left_neighbor -> get_Left_Neighbor() -> get_Node_Location();
    if (left_left_nbr == left_nbr){
        cout << "Second node! My left nbr (1st node) has no angle so F^L = 0!!" << endl;
        return F_left;
    }

    //Now we can deal with the fact that theta_left = 0 or theta_left = theta_equi
    double left_angle = left_neighbor->get_Current_Angle();
    if (left_angle == 0){
        //Has no angle left nbr
        cout << "No angle and No left neighbor! F^L = 0!" << endl;
        return F_left;
    }
    else if (abs(left_angle - pi) < epsilon){
        //takes care of both numerator and denominator being = 0
        cout << "|left_angle - pi| < epsilon! F^L = 0!" << endl;
        return F_left;
    }

    //Now we can calculate F^L for nodes that do have a Left nbr with an angle:
    double left_equi_angle = left_neighbor->get_Equi_Angle();
    double k_bend = left_neighbor->get_K_Bend_Actin();

    constant_left = k_bend*(left_angle - left_equi_angle)/(sqrt(1-pow(cos(left_angle),2)));
    cout << "constant_left = " << constant_left << endl;

    Coord left_vec = left_neighbor -> get_Node_Location() - my_location;
    Coord left_left_vec = left_neighbor -> get_Left_Neighbor() -> get_Node_Location() - left_neighbor -> get_Node_Location();
    double left_length = left_vec.length();
    double left_left_length = left_left_vec.length();

    Coord term_1L = left_left_vec/(left_left_length*left_length);
    Coord term_2L = left_vec*cos(left_angle)/pow(left_length,2); 

	F_left = (term_1L + term_2L) * constant_left;

    cout << "Bending left: " << endl;
    cout << "F_left = (" << term_1L << " + " << term_2L << ") * " << constant_left << " = " << F_left << endl;

    return F_left;

}

Coord Actin_Node::bending_Force_Equation_Right(){
    Coord F_right;
    double constant_right;
    double epsilon = 0.0001;

    //First deal with the fact that the last node has NO angle and NO right nbr (set to be itself)
    Coord right_nbr = right_neighbor -> get_Node_Location();
    if (right_nbr == my_location){
        cout << "Last node! Have NO angle and NO right neighbor so F^R = 0!!" << endl;
        return F_right;
    }

    //The last node has no angle, no right neighbor (set to be itself) and no RR neighbor
    //This is for the 2nd to last node, where F^R = 0 because it's R nbr (last node) has NO angle
    Coord right_right_nbr = right_neighbor -> get_Right_Neighbor() -> get_Node_Location();
    if (right_right_nbr == right_nbr){
        cout << "Second to last node! My right nbr (last node) has no angle so F^R = 0!!" << endl;
        return F_right;
    }

    //Now we can deal with the fact that theta_right = 0 or theta_right = theta_equi
    double right_angle = right_neighbor->get_Current_Angle();
    if (right_angle == 0){
        //Has no angle
        cout << "No angle and No right neighbor! F^R = 0!" << endl;
        return F_right;
    }
    else if (abs(right_angle - pi) < epsilon){
        //takes care of both numerator and denominator being = 0
        cout << "|right_angle - pi| < epsilon! F^R = 0!" << endl;
        return F_right;
    }

    //Now we can calculate F^R for nodes that do have a right nbr with an angle:
    double right_equi_angle = right_neighbor->get_Equi_Angle();
    double k_bend = right_neighbor->get_K_Bend_Actin();

    constant_right = k_bend*(right_angle - right_equi_angle)/(sqrt(1-pow(cos(right_angle),2)));
    cout << "constant right = " << constant_right << endl;

    Coord right_vec = right_neighbor -> get_Node_Location() - my_location;
    Coord right_right_vec = right_neighbor -> get_Right_Neighbor() -> get_Node_Location() - right_neighbor -> get_Node_Location(); 
    double right_length = right_vec.length();
    double right_right_length = right_right_vec.length();

    Coord term_1R = right_right_vec/(right_right_length*right_length);
    Coord term_2R = right_vec*cos(right_angle)/pow(right_length,2);

    F_right = (term_1R + term_2R)*constant_right;

    cout << "Bending right: " << endl;
    cout << "F_right = (" << term_1R << " + " << term_2R << ") * " << constant_right << " = " << F_right << endl;

    return F_right;

}

//***Stochastic force
Coord Actin_Node::calc_Stochastic_Force(){
    Coord F_random;

    //Mean and standard deviation of normal distribution
    double mean = 0.0;
    double stddev = 1.0;

    //generate coord with normally distributed random numbers
    double rand_num1 = get_Random_Number(mean, stddev);
    cout << "rand_num1 = " << rand_num1 << endl;
    double rand_num2 = get_Random_Number(mean, stddev);
    cout << "rand_num2 = " << rand_num2 << endl;

    Coord force = Coord(rand_num1, rand_num2);
    cout << "random coord" << force << endl;

    //now calculate the force equation
    double drag = get_Drag_Coeff();
    double constant = sqrt((2*kB*TEMPERATURE*drag)/(dt));

    F_random = force*constant;
    //F_random = force*sqrt((2*kB*TEMPERATURE*drag)/(dt));

    cout << "F_random = " << force << " * " << constant << " = " << F_random << endl;
    
    return F_random;
}

double Actin_Node::get_Random_Number(double mean, double stddev){
    double rand_num;
    // cout << "mean = " << mean << endl;
    // cout << "stddev = " << stddev << endl;

    rand_num = this->get_My_Filament()->get_Network()->get_Normally_Distributed_Random_Number(mean, stddev);
    
    return rand_num;
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
Myosin_Node::Myosin_Node(Coord location, shared_ptr<Myosin> my_myosin_minifilament) : Node(location){
    this->my_myosin_minifilament = my_myosin_minifilament;
}

Myosin_Node::Myosin_Node(Coord location, shared_ptr<Myosin> my_myosin_minifilament, shared_ptr<Myosin_Node> nbring_pair) : Node(location){
    this->my_myosin_minifilament = my_myosin_minifilament;
    this-> neighboring_pair = nbring_pair;
}

//==============================
//Getters and setters:
//==============================
void Myosin_Node::set_My_Myosin(shared_ptr<Myosin> minifilament){
    this->my_myosin_minifilament = minifilament;
    return;
}

void Myosin_Node::set_Neighboring_Pair(shared_ptr<Myosin_Node> nbh_pair){
    this->neighboring_pair = nbh_pair;
    return;
}

void Myosin_Node::set_K_Linear_Myosin(double k_linear){
    this->k_linear_myosin = k_linear;
    return;
}

void Myosin_Node::set_Myosin_Equi_Len(double equi_len){
    this->myosin_spring_equi_len = equi_len;
    return;
}

//==============================
//Functions:
//==============================
void Myosin_Node::sound_Off_Neighboring_Pair(){
    cout << "Myosin node #: " << get_My_Node_Rank() << endl;
    cout << "Location: " << get_Node_Location() << endl;
    cout << "My neighboring pair: " << get_Neighboring_Pair()->get_Node_Location() << endl;

    return;
}

//==============================
//Destructor:
//==============================
Myosin_Node::~Myosin_Node(){
    cout << "Myosin node destructor called!" << endl;
}

//==========================================================
// End of node.cpp