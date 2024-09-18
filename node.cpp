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

void Actin_Node::set_Connected_Myosin_Node(shared_ptr<Myosin_Node> myo_node){
    this->curr_conn_myo_node = myo_node;
    return;
}

void Actin_Node::set_Conn_Myosin_MiniFilament(shared_ptr<Myosin> minifilament){
    this->curr_conn_minifilament = minifilament;
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
    
    // double crossProduct = left_vec.cross(right_vec);

    // if(crossProduct < 0.0){
    //     theta = 2 * pi - theta;
    // }

    //update the angle and cross product in protected member variables
    my_current_angle = theta;
    cout << "My angle: " << theta << endl;
    //cross_product = crossProduct;

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

void Actin_Node::set_K_Linear_Actomyo_Conn(double k_lin){
    this->k_linear_actomyoConn = k_lin;
    return;
}

void Actin_Node::set_Actomyo_Equi_Len(double equil_len){
    this->actomyo_spring_equi_len = equil_len;
    return;
}

void Actin_Node::set_Myo_Pull_Force(double pull_force){
    this->myo_pull_force = pull_force;
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
    cout << "Calculating total actin force..." << endl;
    cout << "Ti = " << Ti << endl;

    //------Linear force from left and right neighbors------
    //since 1st node doesn't have L nhbr => F^L = 0. Since last node doesn't have R nhbr => F^R = 0.
    //else F^spring = F^L + F^R
    cout << "Linear spring force being calculated..." << endl;
    cout << "==========================================" << endl;
    Coord F_linear = calc_Linear_Force();
    cout << "F_linear = " << F_linear << endl;

    //------Bending force from left and right neighbors------
    cout << "Bending force being calculated..." << endl;
    cout << "==========================================" << endl;
    Coord F_bending = calc_Bending_Force();
    cout << "F_bending = " << F_bending << endl;

    //------Stochastic/random force------
    cout <<"Stochastic/random force being calculated..." << endl;
    cout << "==========================================" << endl;
    Coord F_stoch;
    if(STOCHASTIC_FORCE_ACTIN){
        F_stoch = calc_Stochastic_Force();
        cout << "F_stoch = " << F_stoch << endl;
    }
    else{
        cout << "Stochastic force OFF" << endl;
        cout << "F_stoch = " << F_stoch << endl;
    }
    // Coord F_stoch = calc_Stochastic_Force();
    // cout << "F_stoch = " << F_stoch << endl;

    //------Actomyosin connection force------
    cout <<"Actomyosin connection force being calculated..." << endl;
    cout << "==========================================" << endl;
    Coord F_actomyo_conn;
    if(Ti==0){
        //Coord F_actomyo_conn;
        cout << "No connection! F_actomyo_conn = " << F_actomyo_conn << endl;
    }
    else{
        F_actomyo_conn = calc_Actomyo_Conn_Force();
        cout << "F_actomyo_conn = " << F_actomyo_conn << endl;
    }

    //------Myosin pulling force------
    cout <<"Myosin pulling force due conn being calculated..." << endl;
    cout << "=================================================" << endl;
    Coord F_myo_pull;
    if(Ti==0){
        //Coord F_actomyo_conn;
        cout << "No connection! F_myo_pull = " << F_myo_pull << endl;
    }
    else{
        F_myo_pull = calc_Myosin_Pull_Force();
        cout << "F_myo_pull = " << F_myo_pull << endl;
    }

    //total force actin on node
    new_total_force = F_linear + F_bending + F_stoch + F_actomyo_conn + F_myo_pull;

    cout << "Total force = " << F_linear << " + " << F_bending << " + " << F_stoch << " + " << F_actomyo_conn << " + " << F_myo_pull << " = " << new_total_force << endl;

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

//***Actomyosin connection linear spring force
Coord Actin_Node::calc_Actomyo_Conn_Force(){
    Coord F_actomyo;

    //For double checking-purposes
    cout << "Actin filament # = " << this->get_My_Filament()->get_Filament_Num() << endl;
    cout << "Actin node loc = " << this->get_Node_Location() << endl;

    //Get connected myosin node:
    if(curr_conn_myo_node == NULL){
        cout << "NO ACTOMYO CONNECTION available!" << endl;
        //Will return F_actomyo = (0,0)
    }
    else{
        cout << "Connection with mini-filament = " << curr_conn_minifilament->get_Myosin_Num() << endl;
        cout << "Connected myo node = " << curr_conn_myo_node->get_Node_Location() << endl;

        cout << "Calc linear force from connection with myosin node" << endl;
        Coord F_conn = linear_Actomyo_Spring_Equation(curr_conn_myo_node);
        cout << "F_conn = " << F_conn<< endl;

        F_actomyo = F_conn;
    }

    cout << "F_actomyo = " << F_actomyo << endl;

    return F_actomyo;
}

Coord Actin_Node::linear_Actomyo_Spring_Equation(shared_ptr<Myosin_Node> myo_node){
    if(myo_node == NULL){
        cout << "ERROR: Accessing a NULL pointer. Will abort!!!" << endl;
        exit(1);
    }

    Coord F_connection;

    //compute the difference vector and length
    Coord diff_vec = myo_node->get_Node_Location() - my_location;
    double diff_length = diff_vec.length();
    cout << "diff_vec = " << myo_node->get_Node_Location() << " - " << my_location << " = " << diff_vec << endl;
    cout << "diff_length = " << diff_length << endl;

    //Avoid dividing by 0:
    if(diff_length == 0){
        cout << "Avoiding 0 in the denominator. Returning (0,0)" << endl;
        return Coord(0,0);
    }

    //F_actomyo = k(|node_i - node_j| - l_equi)((node_i - node_j)/|node_i - node_j|)
    cout << "k_linear_actomyoConn = " << k_linear_actomyoConn << endl;
    cout << "actomyo_spring_equi_len = " << actomyo_spring_equi_len << endl;

    F_connection = (diff_vec/diff_length)*k_linear_actomyoConn*(diff_length - this->actomyo_spring_equi_len);
    cout << "Actomyosin connection force: " << F_connection << endl;

    return F_connection;

}

//***Myosin pulling force from conn
Coord Actin_Node::calc_Myosin_Pull_Force(){
    Coord F_pull;

    //For double checking-purposes
    cout << "Actin filament # = " << this->get_My_Filament()->get_Filament_Num() << endl;
    cout << "Actin node loc = " << this->get_Node_Location() << endl;

    //Get connected myosin node:
    if(curr_conn_myo_node == NULL){
        cout << "NO MYO PULLING FORCE because no actomyo connection available!" << endl;
        //Will return F_actomyo = (0,0)
    }
    else{
        cout << "Connection with mini-filament = " << curr_conn_minifilament->get_Myosin_Num() << endl;
        cout << "Connected myo node = " << curr_conn_myo_node->get_Node_Location() << endl;

        cout << "Applying myosin pulling force due to connection..." << endl;
        Coord F_myo = myosin_Pulling_Force(curr_conn_myo_node);
        cout << "F_myo = " << F_myo<< endl;

        F_pull = F_myo;
    }

    cout << "F_pull = " << F_pull << endl;

    return F_pull;
}

Coord Actin_Node::myosin_Pulling_Force(shared_ptr<Myosin_Node> myo_node){
    if(myo_node == NULL){
        cout << "ERROR: Accessing a NULL pointer. Will abort!!!" << endl;
        exit(1);
    }

    Coord F_myosin_force;

    //The pulling force depends on the barbed end, since that determines in which direction the myosin will move.
    //The actin filament slides in the opposite direction
    bool first_isBarbedEnd = this->get_My_Filament()->get_First_Node_Polarity();
    bool last_isBarbedEnd = this->get_My_Filament()->get_Last_Node_Polarity();

    //Get the actin nodes in the filament of this actin
    vector<shared_ptr<Actin_Node>> actins; //will store actin nodes here
    this->get_My_Filament()->get_Actin_Nodes_Vec(actins);

    //For double-checking purposes (comment out once done)
    for(unsigned int k = 0; k < actins.size(); k++){
        Coord actin_location = actins.at(k)->get_Node_Location();
        cout << "actin node = " << actin_location << endl;
    }

    //Get first and last node
    int last = actins.size() - 1;
    shared_ptr<Actin_Node> first_node = actins.at(0);
    shared_ptr<Actin_Node> last_node = actins.at(last);
    Coord first_loc = first_node->get_Node_Location();
    Coord last_loc = last_node->get_Node_Location();

    cout << "first actin node = " << first_loc << endl;
    cout << "last actin node = " << last_loc << endl;
    //cout << "my loc = " << my_location << endl;

    //The barbed end is the FIRST node in the filament: (+) o----o----o----o (-)
        //(*)Filament needs to slide towards the right
        //(*)Myosin to the left to the positive end
    if(first_isBarbedEnd){
        cout << "FIRST NODE is barbed end" << endl;

        if(first_loc == my_location){
            cout << "MADE IT TO BARBED END!!" << endl;

            //first node doesn't have a left neighbor (its itself) so we use right neighbor
            Coord right_nbr = this->get_Right_Neighbor()->get_Node_Location();
            cout << "right neighbor = " << right_nbr << endl;

            //compute the difference vector and length
            Coord diff_vec = right_nbr - my_location;
            double diff_length = diff_vec.length();
            cout << "diff_vec = " << right_nbr << " - " << my_location << " = " << diff_vec << endl;
            cout << "diff_length = " << diff_length << endl;

            //Avoid dividing by 0:
            if(diff_length == 0){
                cout << "Avoiding 0 in the denominator. Returning (0,0)" << endl;
                return Coord(0,0);
            }

            //myosin_force = F^0_myo((node_i - node_j)/|node_i - node_j|)
            cout << "myo_pull_force = " << myo_pull_force << endl;
            
            F_myosin_force = (diff_vec/diff_length)*myo_pull_force;
            
            cout << "F_myosin_force = " << F_myosin_force << endl;

        }
        else{
            //Get left nhbr
            Coord left_nbr = this->get_Left_Neighbor()->get_Node_Location();
            cout << "left neighbor = " << left_nbr << endl;

            //compute the difference vector and length
            Coord diff_vec = my_location - left_nbr;
            double diff_length = diff_vec.length();
            cout << "diff_vec = " << my_location << " - " << left_nbr << " = " << diff_vec << endl;
            cout << "diff_length = " << diff_length << endl;

            //Avoid dividing by 0:
            if(diff_length == 0){
                cout << "Avoiding 0 in the denominator. Returning (0,0)" << endl;
                return Coord(0,0);
            }

            //myosin_force = F^0_myo((node_i - node_j)/|node_i - node_j|)
            cout << "myo_pull_force = " << myo_pull_force << endl;
            
            F_myosin_force = (diff_vec/diff_length)*myo_pull_force;
            
            cout << "F_myosin_force = " << F_myosin_force << endl;
        }

    }
    else if(last_isBarbedEnd){
        //The barbed end is the LAST node in the filament: (-) o----o----o----o (+)
            //(*)Filament needs to slide towards the left
            //(*)Myosin to the right to the positive end   
        cout << "LAST NODE is barbed end" << endl;

        if(last_loc == my_location){
            cout << "MADE IT TO BARBED END!!" << endl;

            //The last node doesn't have a right neighbor (its itself) so we use left neighbor
            Coord left_nbr = this->get_Left_Neighbor()->get_Node_Location();
            cout << "left neighbor = " << left_nbr << endl;

            //compute the difference vector and length
            Coord diff_vec = left_nbr - my_location;
            double diff_length = diff_vec.length();
            cout << "diff_vec = " << left_nbr << " - " << my_location << " = " << diff_vec << endl;
            cout << "diff_length = " << diff_length << endl;

            //Avoid dividing by 0:
            if(diff_length == 0){
                cout << "Avoiding 0 in the denominator. Returning (0,0)" << endl;
                return Coord(0,0);
            }

            //myosin_force = F^0_myo((node_i - node_j)/|node_i - node_j|)
            cout << "myo_pull_force = " << myo_pull_force << endl;
            
            F_myosin_force = (diff_vec/diff_length)*myo_pull_force;
            
            cout << "F_myosin_force = " << F_myosin_force << endl;
        }
        else{
            //get right neighbor:
            Coord right_nbr = this->get_Right_Neighbor()->get_Node_Location();
            cout << "right neighbor = " << right_nbr << endl;

            //compute the difference vector and length
            Coord diff_vec = my_location - right_nbr;
            double diff_length = diff_vec.length();
            cout << "diff_vec = " << my_location << " - " << right_nbr << " = " << diff_vec << endl;
            cout << "diff_length = " << diff_length << endl;

            //Avoid dividing by 0:
            if(diff_length == 0){
                cout << "Avoiding 0 in the denominator. Returning (0,0)" << endl;
                return Coord(0,0);
            }

            //myosin_force = F^0_myo((node_i - node_j)/|node_i - node_j|)
            cout << "myo_pull_force = " << myo_pull_force << endl;
            
            F_myosin_force = (diff_vec/diff_length)*myo_pull_force;
            
            cout << "F_myosin_force = " << F_myosin_force << endl;

        }
    }

    cout << "F_myosin_force = " << F_myosin_force << endl;

    return F_myosin_force;

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

void Myosin_Node::set_If_Connected(bool conn){
    this->connected = conn;
    return;
}

void Myosin_Node::set_Connected_Filament(shared_ptr<Filament> filament){
    this->curr_conn_filament = filament;
    return;
}

void Myosin_Node::set_Connected_Actin_Node(shared_ptr<Actin_Node> node){
    this->curr_conn_actin_node = node;
    return;
}

void Myosin_Node::set_Prev_Connected_Filament(shared_ptr<Filament> prev_filament){
    this->prev_conn_filament = prev_filament;
    return;
}

void Myosin_Node::set_Prev_Connected_Actin_Node(shared_ptr<Actin_Node> prev_node){
    this->prev_conn_actin_node = prev_node;
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

void Myosin_Node::set_K_Linear_Actomyo_Conn(double k_actomyo){
    this->k_linear_actomyo_conn = k_actomyo;
    return;
}

void Myosin_Node::set_Actomyo_Spring_Equi_Length(double equi_len){
    this->actomyo_spring_equi_length = equi_len;
    return;
}

void Myosin_Node::set_Myo_Pulling_Force(double pulling_force){
    this->myo_pulling_force = pulling_force;
    return;
}

//==============================
//Functions:
//==============================
void Myosin_Node::sound_Off_Nyosin_Node_Info(){
    cout << "Parent mini-filament #: " << get_My_Myosin()->get_Myosin_Num() << endl;
    cout << "Node #: " << get_My_Node_Rank() << endl;
    cout << "Location: " << get_Node_Location() << endl;
    cout << "Neighboring pair: " << get_Neighboring_Pair()->get_Node_Location() << endl;
    cout << "Drag coeff: " << get_Drag_Coeff() << endl;
    cout << "Linear spring coeff: " << get_K_Linear_Myosin() << endl;
    cout << "Linear equilibrium length: " << get_Myosin_Equi_Len() << endl;
    cout << "Total force: " << get_Total_Force() << endl;
    cout << "Connected?:" << get_If_Connected() << endl;

    return;
}

void Myosin_Node::sound_Off_Neighboring_Pair(){
    cout << "Myosin node #: " << get_My_Node_Rank() << endl;
    cout << "Location: " << get_Node_Location() << endl;
    cout << "My neighboring pair: " << get_Neighboring_Pair()->get_Node_Location() << endl;

    return;
}

//***Calculate force functions***//
//=================================
//***calculates total force on the myosin node
void Myosin_Node::calculate_Myosin_Forces(int Ti){
    cout << "Calculating total myosin force..." << endl;
    cout << "Ti = " << Ti << endl;

    //------Linear force from neighboring pair------
    cout << "Myosin Linear Spring Force being calculated..." << endl;
    cout << "=================================================" << endl;
    Coord F_linear_myo = calc_Myosin_Linear_Force();
    cout << "F_linear_myo = " << F_linear_myo << endl;

    //------Stochastic/random force------
    cout <<"Myosin Stochastic/random force being calculated..." << endl;
    cout << "=================================================" << endl;
    Coord F_stoch_myo;
    if(STOCHASTIC_FORCE_MYOSIN){
        F_stoch_myo = calc_Myosin_Stochastic_Force();
        cout << "F_stoch_myo = " << F_stoch_myo << endl;
    }
    else{
        cout << "Stochastic force OFF" << endl;
        cout << "F_stoch_myo = " << F_stoch_myo << endl;
    }
    // Coord F_stoch_myo = calc_Myosin_Stochastic_Force();
    // cout << "F_stoch_myo = " << F_stoch_myo << endl;

    //------Actomyosin connection force------
    cout <<"Actomyosin connection force being calculated..." << endl;
    cout << "==========================================" << endl;
    Coord F_actomyo_conn;
    if(Ti==0){
        //Coord F_actomyo_conn;
        cout << "No connection! F_actomyo_conn = " << F_actomyo_conn << endl;
    }
    else{
        F_actomyo_conn = calc_Actomyo_Connection_Force();
        cout << "F_actomyo_conn = " << F_actomyo_conn << endl;
    }

    //------Myosin pulling force------
    cout <<"Myosin pulling force due to conn being calculated..." << endl;
    cout << "=================================================" << endl;
    Coord F_myo_pull;
    if(Ti==0){
        //Coord F_actomyo_conn;
        cout << "No connection! F_myo_pull = " << F_myo_pull << endl;
    }
    else{
        F_myo_pull = calc_Myosin_Pulling_Force();
        cout << "F_myo_pull = " << F_myo_pull << endl;
    }

    //total force myosin on node
    new_total_force = F_linear_myo + F_stoch_myo + F_actomyo_conn + F_myo_pull;

    cout << "Total myosin force = " << F_linear_myo << " + " << F_stoch_myo << " + " << F_actomyo_conn << " + " << F_myo_pull << " = " << new_total_force << endl;

    return;
}

//***Linear spring force
Coord Myosin_Node::calc_Myosin_Linear_Force(){
    Coord F_linear;

    cout << "Calc myosin linear from neighboring pair" << endl;
    Coord F_lin_myo = linear_Myosin_Spring_Equation(neighboring_pair);
    cout << "F_lin_myo = " << F_lin_myo << endl;

    F_linear = F_lin_myo;
    cout << "F_linear_myo = " << F_linear << endl;

    return F_linear;
}

Coord Myosin_Node::linear_Myosin_Spring_Equation(shared_ptr<Myosin_Node> node){
    if(node == NULL){
        cout << "ERROR: Accessing a NULL pointer. Will abort!!!" << endl;
        exit(1);
    }

    Coord F_linear_myo;

    //Compute the vector and length
    Coord diff_vec = node->get_Node_Location() - my_location;
    double diff_length = diff_vec.length();
    cout << "diff_vec_myo = " << diff_vec << endl;
    cout << "diff_length_myo = " << diff_length << endl;

    //To avoid dividing by 0
    if(diff_length == 0){
        cout << "Avoiding 0 in the denominator. Returning (0,0)" << endl;
        return Coord(0,0);
    }

    //F_spring = k(|node_i - node_j| - l_equi)((node_i - node_j)/|node_i - node_j|)
    F_linear_myo = (diff_vec/diff_length)*k_linear_myosin*(diff_length - this->myosin_spring_equi_len);
    cout << "Myosin Linear Force: " << F_linear_myo << endl;

    return F_linear_myo;

}

//***Stochastic force
Coord Myosin_Node::calc_Myosin_Stochastic_Force(){
    Coord F_random_myo;

    //Mean and standard deviation of normal distribution
    double mean = 0.0;
    double stddev = 1.0;

    //Generate coord with normally distributed random numbers
    double myo_rand_num1 = get_Random_Num(mean, stddev);
    cout << "myo_rand_num1 = " << myo_rand_num1 << endl;
    double myo_rand_num2 = get_Random_Num(mean, stddev);
    cout << "myo_rand_num2 = " << myo_rand_num2 << endl;

    Coord rand_force = Coord(myo_rand_num1, myo_rand_num2);
    cout << "myo random coord" << rand_force << endl;

    //Now calculate the random force equation
    double drag = get_Drag_Coeff();
    double constant = sqrt((2*kB*TEMPERATURE*drag)/(dt));

    F_random_myo = rand_force*constant;

    cout << "F_random_myo = " << rand_force << " * " << constant << " = " << F_random_myo << endl;

    return F_random_myo;
}

double Myosin_Node::get_Random_Num(double mean, double stddev){
    double random_number;
    random_number = this->get_My_Myosin()->get_Network()->get_Normally_Distributed_Random_Number(mean, stddev);
    
    return random_number;
}

void Myosin_Node::connect_To_Filament(double radius){
    cout << "connection radius = " << radius << endl;
    cout << "myosin node loc = " << this->get_Node_Location() << endl;

    //Get the neighboring actin filaments
    vector<shared_ptr<Filament>> nbhrActinFilaments;
    nbhrActinFilaments = this -> get_My_Myosin() -> get_Possible_Connections();
    int num_Possible_Fils = nbhrActinFilaments.size();

    //Get if myosin node is connected
    bool isConnected = this->get_If_Connected();

    if(nbhrActinFilaments.empty()){
        cout << "NO neighboring actin filaments. NO actomyo connection formed!!" << endl;
    }
    else{
        //for double checking purposes:
        cout << "Neighboring actin filaments:" << endl;
        for(unsigned int i = 0; i < nbhrActinFilaments.size(); i++){
            cout << nbhrActinFilaments.at(i)->get_Filament_Num() << endl;
        }

        //If a connection exists, see if nbr actin node in the + direction is in range and update connection. If not keep connection. Also check if made it to barbed end
        if(isConnected){
            cout << "CURRENTLY CONNECTED. Myosin head looking to update connection..." << endl;

            //Get myosin's connection info:
            shared_ptr<Filament> connected_Filament = this->get_Connected_Filament();
            shared_ptr<Actin_Node> connected_Actin_Node = this->get_Connected_Actin_Node();

            int conn_Actin_Fil = connected_Filament->get_Filament_Num();
            Coord conn_Actin_Node = connected_Actin_Node->get_Node_Location();

            cout << "CURR connected filament = " << conn_Actin_Fil << endl;
            cout << "CURR connected actin node = " << conn_Actin_Node << endl;

            //Get the connected filament's polarity
            bool actin_First_isBarbed = connected_Filament->get_First_Node_Polarity();
            bool actin_Last_isBarbed = connected_Filament->get_Last_Node_Polarity();

            cout << "FIRST is barbed end?: " << actin_First_isBarbed << endl;
            cout << "LAST is barbed end?: " << actin_Last_isBarbed << endl;

            //Get the connected filament's nodes
            vector<shared_ptr<Actin_Node>> actins; //will store actin nodes here
            connected_Filament->get_Actin_Nodes_Vec(actins);

            //For double-checking purposes (comment out once done)
            for(unsigned int k = 0; k < actins.size(); k++){
                Coord actin_location = actins.at(k)->get_Node_Location();
                cout << "actin node = " << actin_location << endl;
            }

            //Get location of the first and last actin node of connected filament
            int last = actins.size() - 1;
            Coord first_Actin_Node = actins.at(0)->get_Node_Location();
            Coord last_Actin_Node = actins.at(last)->get_Node_Location();

            cout << "FIRST actin node = " << first_Actin_Node << endl;
            cout << "LAST actin node = " << last_Actin_Node << endl;

            //FIRST check if made it to barbed end. If reached (+) end
                //1) disconnect i.e. connected = 0, clear current connected node and filament for both actin and myosin
                //2) possibly remove filament from pool of possible filaments

            if(actin_First_isBarbed == true && (conn_Actin_Node == first_Actin_Node)){
            //if(actin_First_isBarbed == true){
                cout << "MADE IT TO THE BARBED END (FIRST node is barbed)!!" << endl;
                cout << "Disconnecting..." << endl;

                //First store old connections:
                this->set_Prev_Connected_Filament(this->get_Connected_Filament());
                this->set_Prev_Connected_Actin_Node(this->get_Connected_Actin_Node());

                cout << "PREV connected filament = " << this->get_Prev_Connected_Filament()-> get_Filament_Num() << endl;
                cout << "PREV connected actin node = " << this->get_Prev_Connected_Actin_Node()->get_Node_Location() << endl;

                //First clear actin node connections
                this->get_Connected_Actin_Node()->set_Conn_Myosin_MiniFilament(nullptr);
                this->get_Connected_Actin_Node()->set_Connected_Myosin_Node(nullptr);

                cout << "actin node loc = " << conn_Actin_Node << endl;
                cout << "verifying loc = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;
                cout << "actin connected to minifilament = " << (this->get_Connected_Actin_Node()->get_Conn_Myosin_MiniFilament() ? "Not null" : "null") << endl;
                cout << "actin connected to node = " << (this->get_Connected_Actin_Node()->get_Conn_Myosin_Node() ? "Not null" : "null") << endl;
                
                //Disconnect and clear myosin node connections:
                this->set_If_Connected(false); //connected = 0
                this->set_Connected_Filament(nullptr);
                this->set_Connected_Actin_Node(nullptr);

                cout << "MYOSIN node connected? = " << this->get_If_Connected() << endl;
                cout << "connected filament = " << (this->get_Connected_Filament() ? "Not null" : "null") << endl;
                cout << "connected actin node = " << (this->get_Connected_Actin_Node()? "Not null" : "null") << endl;

                cout << "DISCONNECTED!" << endl;
            }
            else if(actin_Last_isBarbed == true && (conn_Actin_Node == last_Actin_Node)){
            //else if(actin_Last_isBarbed == true){
                cout << "MADE IT TO THE BARBED END (LAST node is barbed)!!" << endl;
                cout << "Disconnecting..." << endl;

                //First store old connections:
                this->set_Prev_Connected_Filament(this->get_Connected_Filament());
                this->set_Prev_Connected_Actin_Node(this->get_Connected_Actin_Node());

                cout << "PREV connected filament = " << this->get_Prev_Connected_Filament()-> get_Filament_Num() << endl;
                cout << "PREV connected actin node = " << this->get_Prev_Connected_Actin_Node()->get_Node_Location() << endl;

                //First clear actin node connections
                this->get_Connected_Actin_Node()->set_Conn_Myosin_MiniFilament(nullptr);
                this->get_Connected_Actin_Node()->set_Connected_Myosin_Node(nullptr);

                cout << "actin node loc = " << conn_Actin_Node << endl;
                cout << "verifying loc = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;
                cout << "actin connected to minifilament = " << (this->get_Connected_Actin_Node()->get_Conn_Myosin_MiniFilament() ? "Not null" : "null") << endl;
                cout << "actin connected to node = " << (this->get_Connected_Actin_Node()->get_Conn_Myosin_Node() ? "Not null" : "null") << endl;
                
                //Disconnect and clear myosin node connections:
                this->set_If_Connected(false); //connected = 0
                this->set_Connected_Filament(nullptr);
                this->set_Connected_Actin_Node(nullptr);

                cout << "MYOSIN node connected? = " << this->get_If_Connected() << endl;
                cout << "connected filament = " << (this->get_Connected_Filament() ? "Not null" : "null") << endl;
                cout << "connected actin node = " << (this->get_Connected_Actin_Node()? "Not null" : "null") << endl;

                cout << "DISCONNECTED!" << endl;
            }
            else{
                //NEXT, see if nbr actin node in the + direction is in range and update connection. If not keep connection.
                cout << "SEARCHING if nbring actin node in the + direction is in range..." << endl;

                //If barbed end is FIRST node (move in Left direction)
                if(actin_First_isBarbed){
                    cout << "First actin node is + end..." << endl;

                    //Get current actin connected node's LEFT nhbr
                    shared_ptr<Actin_Node> left_Actin_Nbr = connected_Actin_Node->get_Left_Neighbor();
                    Coord left_Actin_Loc = left_Actin_Nbr->get_Node_Location();

                    cout << "LEFT actin node (possible conn) = " << left_Actin_Loc << endl;
                
                    //Calculate the distance between the myosin node and NEW possible conn
                    Coord diff_vec = left_Actin_Loc - my_location;
                    cout << "diff vec = " << left_Actin_Loc << " - " << my_location << " = " << diff_vec << endl;

                    double distance = diff_vec.length();
                    cout << "distance = " << distance << endl;

                    //double r_capture = 0.35; //for testing while writing code 

                    //Check is the distance is less than or equal to connection_radius. If yes, form connection
                    //if(distance <= r_capture){ //for testing while writing code
                    if(distance <= radius){
                        cout << "NEW CONNECTION found towards (+) end.." << endl;

                        //Store previous connections:
                        this->set_Prev_Connected_Actin_Node(this->get_Connected_Actin_Node());
                        this->set_Prev_Connected_Filament(this->get_Connected_Filament());
                        cout << "PREV connected filament = " << this->get_Prev_Connected_Filament()-> get_Filament_Num() << endl;
                        cout << "PREV connected actin node = " << this->get_Prev_Connected_Actin_Node()->get_Node_Location() << endl;

                        //RELEASE the CURRENT connected node
                        this->get_Connected_Actin_Node()->set_Conn_Myosin_MiniFilament(nullptr);
                        this->get_Connected_Actin_Node()->set_Connected_Myosin_Node(nullptr);

                        cout << "actin node loc = " << conn_Actin_Node << endl;
                        cout << "verifying loc = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;
                        cout << "actin connected to minifilament = " << (this->get_Connected_Actin_Node()->get_Conn_Myosin_MiniFilament() ? "Not null" : "null") << endl;
                        cout << "actin connected to node = " << (this->get_Connected_Actin_Node()->get_Conn_Myosin_Node() ? "Not null" : "null") << endl;
                
                        cout << "UPDATING CONNECTIONS..." << endl;

                        //Update new connections:
                        //MYOSIN:
                        this->set_Connected_Actin_Node(left_Actin_Nbr);
                        this->set_Connected_Filament(connected_Filament);

                        cout << "myosin connected = " << this->get_If_Connected() << endl;
                        cout << "connected filament = " << this->get_Connected_Filament()->get_Filament_Num() << endl;
                        cout << "NEW connected actin node = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;

                        //ACTIN:
                        left_Actin_Nbr->set_Connected_Myosin_Node(this->shared_from_this());
                        left_Actin_Nbr->set_Conn_Myosin_MiniFilament(this->get_My_Myosin());

                        cout << "NEW actin node loc = " << left_Actin_Nbr->get_Node_Location() << endl;
                        cout << "NEW actin connected to node = " << left_Actin_Nbr->get_Conn_Myosin_Node()->get_Node_Location() << endl;
                        cout << "NEW actin connected to minifilament = " << left_Actin_Nbr->get_Conn_Myosin_MiniFilament()->get_Myosin_Num() << endl;

                        cout << "NEW CONNECTION MADE!" << endl;
                            
                    }
                    else{
                        // cout << "KEEPING prev connections..." << endl;

                        // cout << "myosin connected = " << this->get_If_Connected() << endl;
                        // cout << "SAME connected filament = " << this->get_Connected_Filament()->get_Filament_Num() << endl;
                        // cout << "SAME connected actin node = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;

                        // cout << "actin node loc = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;
                        // cout << "actin connected to minifilament = " << this->get_Connected_Actin_Node()->get_Conn_Myosin_MiniFilament()->get_Myosin_Num() << endl;
                        // cout << "actin connected to node = " << this->get_Connected_Actin_Node()->get_Conn_Myosin_Node()->get_Node_Location() << endl;

                        //Check current connection
                        cout << "CHECKING if prev connection still in range..." << endl;
                        cout << "myosin connected = " << this->get_If_Connected() << endl;
                        cout << "CURRENT connected filament = " << this->get_Connected_Filament()->get_Filament_Num() << endl;
                        cout << "CURRENT connected actin node = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;

                        Coord current_conn_node = this->get_Connected_Actin_Node()->get_Node_Location();

                        //Calculate the distance between the myosin node and NEW possible conn
                        Coord diff_vec = current_conn_node - my_location;
                        cout << "diff vec = " << current_conn_node << " - " << my_location << " = " << diff_vec << endl;

                        double distance = diff_vec.length();
                        cout << "distance = " << distance << endl;  

                        //If still in range, KEEP same connection
                        if(distance <= radius){
                            cout << "KEEPING prev connections..." << endl;

                            cout << "myosin connected = " << this->get_If_Connected() << endl;
                            cout << "SAME connected filament = " << this->get_Connected_Filament()->get_Filament_Num() << endl;
                            cout << "SAME connected actin node = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;

                            cout << "actin node loc = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;
                            // cout << "actin connected to minifilament = " << this->get_Connected_Actin_Node()->get_Conn_Myosin_MiniFilament()->get_Myosin_Num() << endl;
                            // cout << "actin connected to node = " << this->get_Connected_Actin_Node()->get_Conn_Myosin_Node()->get_Node_Location() << endl;
                        }
                        else{
                            //If no longer in range, DISCONNECT
                            cout << "PREVIOUS connection NO LONGER in range..." << endl;
                            cout << "Disconnecting..." << endl;

                            //First store old connections:
                            this->set_Prev_Connected_Filament(this->get_Connected_Filament());
                            this->set_Prev_Connected_Actin_Node(this->get_Connected_Actin_Node());

                            cout << "PREV connected filament = " << this->get_Prev_Connected_Filament()-> get_Filament_Num() << endl;
                            cout << "PREV connected actin node = " << this->get_Prev_Connected_Actin_Node()->get_Node_Location() << endl;

                            //First clear actin node connections
                            this->get_Connected_Actin_Node()->set_Conn_Myosin_MiniFilament(nullptr);
                            this->get_Connected_Actin_Node()->set_Connected_Myosin_Node(nullptr);

                            cout << "actin node loc = " << conn_Actin_Node << endl;
                            cout << "verifying loc = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;
                            cout << "actin connected to minifilament = " << (this->get_Connected_Actin_Node()->get_Conn_Myosin_MiniFilament() ? "Not null" : "null") << endl;
                            cout << "actin connected to node = " << (this->get_Connected_Actin_Node()->get_Conn_Myosin_Node() ? "Not null" : "null") << endl;
                            
                            //Disconnect and clear myosin node connections:
                            this->set_If_Connected(false); //connected = 0
                            this->set_Connected_Filament(nullptr);
                            this->set_Connected_Actin_Node(nullptr);

                            cout << "MYOSIN node connected? = " << this->get_If_Connected() << endl;
                            cout << "connected filament = " << (this->get_Connected_Filament() ? "Not null" : "null") << endl;
                            cout << "connected actin node = " << (this->get_Connected_Actin_Node()? "Not null" : "null") << endl;

                            cout << "DISCONNECTED!" << endl;

                        }

                    }

                }
                else if(actin_Last_isBarbed){
                    //If barbed end is LAST node move in Right direction
                    cout << "Last actin node is + end..." << endl;

                    //Get current actin connected node's RIGHT nhbr
                    shared_ptr<Actin_Node> right_Actin_Nbr = connected_Actin_Node->get_Right_Neighbor();
                    Coord right_Actin_Loc = right_Actin_Nbr->get_Node_Location();

                    cout << "RIGHT actin node (possible conn) = " << right_Actin_Loc << endl;

                    //Calculate the distance between the myosin node and NEW possible conn
                    Coord diff_vec = right_Actin_Loc - my_location;
                    cout << "diff vec = " << right_Actin_Loc << " - " << my_location << " = " << diff_vec << endl;

                    double distance = diff_vec.length();
                    cout << "distance = " << distance << endl;

                    //double r_capture = 0.35; //for testing while writing code 

                    //Check is the distance is less than or equal to connection_radius. If yes, form connection
                    //if(distance <= r_capture){ //for testing while writing code
                    if(distance <= radius){
                        cout << "NEW CONNECTION found towards (+) end.." << endl;

                        //Store previous connections:
                        this->set_Prev_Connected_Actin_Node(this->get_Connected_Actin_Node());
                        this->set_Prev_Connected_Filament(this->get_Connected_Filament());
                        cout << "PREV connected filament = " << this->get_Prev_Connected_Filament()-> get_Filament_Num() << endl;
                        cout << "PREV connected actin node = " << this->get_Prev_Connected_Actin_Node()->get_Node_Location() << endl;

                        //RELEASE the CURRENT connected node
                        this->get_Connected_Actin_Node()->set_Conn_Myosin_MiniFilament(nullptr);
                        this->get_Connected_Actin_Node()->set_Connected_Myosin_Node(nullptr);

                        cout << "actin node loc = " << conn_Actin_Node << endl;
                        cout << "verifying loc = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;
                        cout << "actin connected to minifilament = " << (this->get_Connected_Actin_Node()->get_Conn_Myosin_MiniFilament() ? "Not null" : "null") << endl;
                        cout << "actin connected to node = " << (this->get_Connected_Actin_Node()->get_Conn_Myosin_Node() ? "Not null" : "null") << endl;
                    
                        cout << "UPDATING CONNECTIONS..." << endl;

                        //Update new connections:
                        //MYOSIN:
                        this->set_Connected_Actin_Node(right_Actin_Nbr);
                        this->set_Connected_Filament(connected_Filament);

                        cout << "myosin connected = " << this->get_If_Connected() << endl;
                        cout << "connected filament = " << this->get_Connected_Filament()->get_Filament_Num() << endl;
                        cout << "NEW connected actin node = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;

                        //ACTIN:
                        right_Actin_Nbr->set_Connected_Myosin_Node(this->shared_from_this());
                        right_Actin_Nbr->set_Conn_Myosin_MiniFilament(this->get_My_Myosin());

                        cout << "NEW actin node loc = " << right_Actin_Nbr->get_Node_Location() << endl;
                        cout << "NEW actin connected to node = " << right_Actin_Nbr->get_Conn_Myosin_Node()->get_Node_Location() << endl;
                        cout << "NEW actin connected to minifilament = " << right_Actin_Nbr->get_Conn_Myosin_MiniFilament()->get_Myosin_Num() << endl;

                        cout << "NEW CONNECTION MADE!" << endl;
                    }
                    else{
                        // cout << "KEEPING prev connections..." << endl;

                        // cout << "myosin connected = " << this->get_If_Connected() << endl;
                        // cout << "SAME connected filament = " << this->get_Connected_Filament()->get_Filament_Num() << endl;
                        // cout << "SAME connected actin node = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;

                        // cout << "actin node loc = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;
                        // cout << "actin connected to minifilament = " << this->get_Connected_Actin_Node()->get_Conn_Myosin_MiniFilament()->get_Myosin_Num() << endl;
                        // cout << "actin connected to node = " << this->get_Connected_Actin_Node()->get_Conn_Myosin_Node()->get_Node_Location() << endl;

                        //Check current connection
                        cout << "CHECKING if prev connection still in range..." << endl;
                        cout << "myosin connected = " << this->get_If_Connected() << endl;
                        cout << "CURRENT connected filament = " << this->get_Connected_Filament()->get_Filament_Num() << endl;
                        cout << "CURRENT connected actin node = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;

                        Coord current_conn_node = this->get_Connected_Actin_Node()->get_Node_Location();

                        //Calculate the distance between the myosin node and NEW possible conn
                        Coord diff_vec = current_conn_node - my_location;
                        cout << "diff vec = " << current_conn_node << " - " << my_location << " = " << diff_vec << endl;

                        double distance = diff_vec.length();
                        cout << "distance = " << distance << endl;  

                        //If still in range, KEEP same connection
                        if(distance <= radius){
                            cout << "KEEPING prev connections..." << endl;

                            cout << "myosin connected = " << this->get_If_Connected() << endl;
                            cout << "SAME connected filament = " << this->get_Connected_Filament()->get_Filament_Num() << endl;
                            cout << "SAME connected actin node = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;

                            cout << "actin node loc = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;
                            // cout << "actin connected to minifilament = " << this->get_Connected_Actin_Node()->get_Conn_Myosin_MiniFilament()->get_Myosin_Num() << endl;
                            // cout << "actin connected to node = " << this->get_Connected_Actin_Node()->get_Conn_Myosin_Node()->get_Node_Location() << endl;
                        }
                        else{
                            //If no longer in range, DISCONNECT
                            cout << "PREVIOUS connection NO LONGER in range..." << endl;
                            cout << "Disconnecting..." << endl;

                            //First store old connections:
                            this->set_Prev_Connected_Filament(this->get_Connected_Filament());
                            this->set_Prev_Connected_Actin_Node(this->get_Connected_Actin_Node());

                            cout << "PREV connected filament = " << this->get_Prev_Connected_Filament()-> get_Filament_Num() << endl;
                            cout << "PREV connected actin node = " << this->get_Prev_Connected_Actin_Node()->get_Node_Location() << endl;

                            //First clear actin node connections
                            this->get_Connected_Actin_Node()->set_Conn_Myosin_MiniFilament(nullptr);
                            this->get_Connected_Actin_Node()->set_Connected_Myosin_Node(nullptr);

                            cout << "actin node loc = " << conn_Actin_Node << endl;
                            cout << "verifying loc = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;
                            cout << "actin connected to minifilament = " << (this->get_Connected_Actin_Node()->get_Conn_Myosin_MiniFilament() ? "Not null" : "null") << endl;
                            cout << "actin connected to node = " << (this->get_Connected_Actin_Node()->get_Conn_Myosin_Node() ? "Not null" : "null") << endl;
                            
                            //Disconnect and clear myosin node connections:
                            this->set_If_Connected(false); //connected = 0
                            this->set_Connected_Filament(nullptr);
                            this->set_Connected_Actin_Node(nullptr);

                            cout << "MYOSIN node connected? = " << this->get_If_Connected() << endl;
                            cout << "connected filament = " << (this->get_Connected_Filament() ? "Not null" : "null") << endl;
                            cout << "connected actin node = " << (this->get_Connected_Actin_Node()? "Not null" : "null") << endl;

                            cout << "DISCONNECTED!" << endl;

                        }
                    }

                }
 

            }
            
        }
        else{
            //If no connection has been formed (connected = 0), search for a connection
            cout << "NOT CONNECTED. Myosin head searching for a connection..." << endl;

            //Get nhbring pairs connected filament. Head#1 and head#2 CANNOT bind to the same actin filament
            shared_ptr<Filament> pairs_Conn_Filament = this->get_Neighboring_Pair()->get_Connected_Filament();
            cout << "Neighboring pairs connected filament = " << pairs_Conn_Filament<< endl;
            
            //Loop through every actin node of every filament to find connection
            for(unsigned int i = 0; i < nbhrActinFilaments.size(); i++){
                int fil_num = nbhrActinFilaments.at(i)->get_Filament_Num();
                cout << "Actin filament #:" << fil_num << endl;

                //get the actin nodes of nhbring filament
                vector<shared_ptr<Actin_Node>> actins; //will store actin nodes here
                nbhrActinFilaments.at(i)->get_Actin_Nodes_Vec(actins);
                //shared_ptr<Actin_Node> actin_fil_node = NULL;

                if(nbhrActinFilaments.at(i) != pairs_Conn_Filament){
                    //Proceed to find connection
                    cout << "Searching through available neighboring actin filaments..." << endl;

                    //Go through every node in the actin filament
                    for(unsigned int k = 0; k < actins.size(); k++){
                        //actin_fil_node = actins.at(k);
                        Coord actin_location = actins.at(k)->get_Node_Location();
                        cout << "actin node = " << actin_location << endl;

                        //Calculate the distance between the actin nodes and myosin node:
                        Coord diff_vec = actin_location - my_location;
                        cout << "diff vec = " << actin_location << " - " << my_location << " = " << diff_vec << endl;

                        double distance = diff_vec.length();
                        cout << "distance = " << distance << endl;

                        //Check is the distance is less than or equal to connection_radius. If yes, form connection
                        if(distance <= radius){
                            cout << "CONNECTION FOUND.." << endl;

                            //set myosin connections
                            this->set_If_Connected(true); //connected = 1
                            this->set_Connected_Filament(nbhrActinFilaments.at(i));
                            this->set_Connected_Actin_Node(actins.at(k));

                            cout << "connected = " << this->get_If_Connected() << endl;
                            cout << "connected filament = " << this->get_Connected_Filament()->get_Filament_Num() << endl;
                            cout << "connected actin node = " << this->get_Connected_Actin_Node()->get_Node_Location() << endl;

                            //set actin connections
                            actins.at(k)->set_Connected_Myosin_Node(this->shared_from_this());
                            actins.at(k)->set_Conn_Myosin_MiniFilament(this->get_My_Myosin());

                            cout << "actin node loc = " << actins.at(k)->get_Node_Location() << endl;
                            cout << "actin connected to node = " << actins.at(k)->get_Conn_Myosin_Node()->get_Node_Location() << endl;
                            cout << "actin connected to minifilament = " << actins.at(k)->get_Conn_Myosin_MiniFilament()->get_Myosin_Num() << endl;

                            cout << "CONNECTION MADE!" << endl;
                            break; //move on to the next actin filament
                        }
                    } 

                    //if a connection was made, do not search through the rest of filaments
                    if(this->get_If_Connected()){
                        break;
                    }
                    
                }
                else if(num_Possible_Fils==1 && nbhrActinFilaments.at(i) == pairs_Conn_Filament){
                    cout << "No connection available! Pairing head is already connected to only possible neghboring filament" << endl;
                }
            }

            bool connctd = this->get_If_Connected();
            if(connctd == false){
                cout << "NO CONNECTION FOUND!!" << endl;
            }

        }
 
    }
    return;
}

//***Actomyosin connection linear spring force
Coord Myosin_Node::calc_Actomyo_Connection_Force(){
    Coord F_actomyo;

    //For double checking-purposes
    cout << "Myosin mini-filament # = " << this->get_My_Myosin()->get_Myosin_Num() << endl;
    cout << "Myosin node loc = " << this->get_Node_Location() << endl;

    if(connected == false || curr_conn_actin_node == NULL){
        cout << "NO ACTOMYO CONNECTION available!" << endl;
        //Will return F_actomyo = (0,0)
    }
    else{
        cout << "Connection with filament = " << curr_conn_filament->get_Filament_Num() << endl;
        cout << "Connected actin node = " << curr_conn_actin_node->get_Node_Location() << endl;

        cout << "Calc linear force from connection with actin node" << endl;
        Coord F_conn = linear_Spring_Actomyo_Equation(curr_conn_actin_node);
        cout << "F_conn = " << F_conn<< endl;

        F_actomyo = F_conn;
    }

    cout << "F_actomyo = " << F_actomyo << endl;

    return F_actomyo;
}

Coord Myosin_Node::linear_Spring_Actomyo_Equation(shared_ptr<Actin_Node> actin_node){
     if(actin_node == NULL){
        cout << "ERROR: Accessing a NULL pointer. Will abort!!!" << endl;
        exit(1);
    }

    Coord F_connection;

    //compute the difference vector and length
    Coord diff_vec = actin_node->get_Node_Location() - my_location;
    double diff_length = diff_vec.length();
    cout << "diff_vec = " << actin_node->get_Node_Location() << " - " << my_location << " = " << diff_vec << endl;
    cout << "diff_length = " << diff_length << endl;

    //Avoid dividing by 0:
    if(diff_length == 0){
        cout << "Avoiding 0 in the denominator. Returning (0,0)" << endl;
        return Coord(0,0);
    }

    //F_actomyo = k(|node_i - node_j| - l_equi)((node_i - node_j)/|node_i - node_j|)
    cout << "k_linear_actomyo_conn = " << k_linear_actomyo_conn << endl;
    cout << "actomyo_spring_equi_length = " << actomyo_spring_equi_length << endl;

    F_connection = (diff_vec/diff_length)*k_linear_actomyo_conn*(diff_length - this->actomyo_spring_equi_length);
    cout << "Actomyosin connection force: " << F_connection << endl;

    return F_connection;

}

//***Myosin pulling force from connection
Coord Myosin_Node::calc_Myosin_Pulling_Force(){
    Coord F_pulling;

    //For double checking-purposes
    cout << "Myosin mini-filament # = " << this->get_My_Myosin()->get_Myosin_Num() << endl;
    cout << "Myosin node loc = " << this->get_Node_Location() << endl;

    if(connected == false || curr_conn_actin_node == NULL){
        cout << "NO ACTOMYO CONNECTION available!" << endl;
        //Will return F_actomyo = (0,0)
    }
    else{
        cout << "Connection with filament = " << curr_conn_filament->get_Filament_Num() << endl;
        cout << "Connected actin node = " << curr_conn_actin_node->get_Node_Location() << endl;

        cout << "Applying myosin pulling force due to connection..." << endl;
        Coord F_myo = myosin_Pulling_Force(curr_conn_filament, curr_conn_actin_node);
        cout << "F_myo = " << F_myo<< endl;

        F_pulling = F_myo;
    }

    cout << "F_pulling = " << F_pulling << endl;

    return F_pulling;
}

Coord Myosin_Node::myosin_Pulling_Force(shared_ptr<Filament> conn_filament, shared_ptr<Actin_Node> actin_node){
    if(actin_node == NULL){
        cout << "ERROR: Accessing a NULL pointer. Will abort!!!" << endl;
        exit(1);
    }

    Coord F_myosin_pull;

    // cout << "Connected filament = " << conn_filament->get_Filament_Num() << endl;
    // cout << "Connected node = " << actin_node->get_Node_Location() << endl;
    
    //Get the actin nodes in the filament of this actin
    vector<shared_ptr<Actin_Node>> actins; //will store actin nodes here
    conn_filament->get_Actin_Nodes_Vec(actins);

    //For double-checking purposes (comment out once done)
    for(unsigned int k = 0; k < actins.size(); k++){
        Coord actin_location = actins.at(k)->get_Node_Location();
        cout << "actin node = " << actin_location << endl;
    }

    //Get the connected filaments polarity:
    //The pulling force depends on the barbed end, since that determines in which direction the myosin will move.
    //The actin filament slides in the opposite direction
    bool first_isBarbedEnd = conn_filament->get_First_Node_Polarity();
    bool last_isBarbedEnd = conn_filament->get_Last_Node_Polarity();

    //Get first and last node
    int last = actins.size() - 1;
    shared_ptr<Actin_Node> first_node = actins.at(0);
    shared_ptr<Actin_Node> last_node = actins.at(last);
    Coord first_loc = first_node->get_Node_Location();
    Coord last_loc = last_node->get_Node_Location();

    cout << "first actin node = " << first_loc << endl;
    cout << "last actin node = " << last_loc << endl;

    //Get the connected filament's nodes
    Coord curr_actin_node = actin_node->get_Node_Location();
    cout << "curr_actin_node = " << curr_actin_node << endl;

    //The barbed end is the FIRST node in the filament: (+) o----o----o----o (-)
        //(*)Filament needs to slide towards the right
        //(*)Myosin to the left to the positive end
    if(first_isBarbedEnd){
        cout << "FIRST ACTIN NODE is barbed end" << endl;

        if(first_loc == curr_actin_node){
            cout << "MADE IT TO BARBED END!!" << endl;

            //First actin node doesn't have a left neighbor (its itself) so we use right neighbor
            Coord right_actin_node = actin_node->get_Right_Neighbor()->get_Node_Location();
            cout << "right actin neighbor = " << right_actin_node << endl;

            //compute the difference vector and length
            Coord diff_vec = curr_actin_node - right_actin_node;
            double diff_length = diff_vec.length();
            cout << "diff_vec = " << right_actin_node << " - " << curr_actin_node << " = " << diff_vec << endl;
            cout << "diff_length = " << diff_length << endl;

            //Avoid dividing by 0:
            if(diff_length == 0){
                cout << "Avoiding 0 in the denominator. Returning (0,0)" << endl;
                return Coord(0,0);
            }

            //myosin_force = F^0_myo((node_i - node_j)/|node_i - node_j|)
            cout << "myo_pulling_force = " << myo_pulling_force << endl;
            
            F_myosin_pull = (diff_vec/diff_length)*myo_pulling_force;
            
            cout << "F_myosin_pull = " << F_myosin_pull << endl;

        }
        else{
            //Get the left neighbor
            Coord left_actin_node = actin_node->get_Left_Neighbor()->get_Node_Location();
            cout << "left actin neighbor = " << left_actin_node << endl;

            //compute the difference vector and length
            Coord diff_vec = left_actin_node - curr_actin_node;
            double diff_length = diff_vec.length();
            cout << "diff_vec = " << left_actin_node << " - " << curr_actin_node << " = " << diff_vec << endl;
            cout << "diff_length = " << diff_length << endl;

            //Avoid dividing by 0:
            if(diff_length == 0){
                cout << "Avoiding 0 in the denominator. Returning (0,0)" << endl;
                return Coord(0,0);
            }

            //myosin_force = F^0_myo((node_i - node_j)/|node_i - node_j|)
            cout << "myo_pulling_force = " << myo_pulling_force << endl;
            
            F_myosin_pull = (diff_vec/diff_length)*myo_pulling_force;
            
            cout << "F_myosin_pull = " << F_myosin_pull << endl;

        }
    }
    else if(last_isBarbedEnd){
        //The barbed end is the LAST node in the filament: (-) o----o----o----o (+)
            //(*)Filament needs to slide towards the left
            //(*)Myosin to the right to the positive end   
        cout << "LAST ACTIN NODE is barbed end" << endl;

        if(last_loc == curr_actin_node){
            cout << "MADE IT TO BARBED END!!" << endl;

            //The last node doesn't have a right neighbor (its itself) so we use left neighbor
            Coord left_actin_node = actin_node->get_Left_Neighbor()->get_Node_Location();
            cout << "left actin neighbor = " << left_actin_node << endl;

            //compute the difference vector and length
            Coord diff_vec = curr_actin_node - left_actin_node;
            double diff_length = diff_vec.length();
            cout << "diff_vec = " << curr_actin_node << " - " << left_actin_node << " = " << diff_vec << endl;
            cout << "diff_length = " << diff_length << endl;

            //Avoid dividing by 0:
            if(diff_length == 0){
                cout << "Avoiding 0 in the denominator. Returning (0,0)" << endl;
                return Coord(0,0);
            }

            //myosin_force = F^0_myo((node_i - node_j)/|node_i - node_j|)
            cout << "myo_pulling_force = " << myo_pulling_force << endl;
            
            F_myosin_pull = (diff_vec/diff_length)*myo_pulling_force;
            
            cout << "F_myosin_pull = " << F_myosin_pull << endl;

        }
        else{
            //Get right neighbor
            Coord right_actin_node = actin_node->get_Right_Neighbor()->get_Node_Location();
            cout << "right actin neighbor = " << right_actin_node << endl;

            //compute the difference vector and length
            Coord diff_vec = right_actin_node - curr_actin_node;
            double diff_length = diff_vec.length();
            cout << "diff_vec = " << right_actin_node << " - " << curr_actin_node << " = " << diff_vec << endl;
            cout << "diff_length = " << diff_length << endl;

            //Avoid dividing by 0:
            if(diff_length == 0){
                cout << "Avoiding 0 in the denominator. Returning (0,0)" << endl;
                return Coord(0,0);
            }

            //myosin_force = F^0_myo((node_i - node_j)/|node_i - node_j|)
            cout << "myo_pulling_force = " << myo_pulling_force << endl;
            
            F_myosin_pull = (diff_vec/diff_length)*myo_pulling_force;
            
            cout << "F_myosin_pull = " << F_myosin_pull << endl;
        }
    }
 
    cout << "F_myosin_pull = " << F_myosin_pull << endl;

    return F_myosin_pull;
}

//==============================
//Destructor:
//==============================
Myosin_Node::~Myosin_Node(){
    cout << "Myosin node destructor called!" << endl;
}

//==========================================================
// End of node.cpp