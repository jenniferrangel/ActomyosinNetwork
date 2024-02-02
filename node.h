//node.h
//=====================
// Include Guards
#ifndef _NODE_H_INCLUDED_  //if node.h hasn't been included yet
#define _NODE_H_INCLUDED_  //    define it so the compiler knows 
//=====================
// Forward Declarations 
class Actin_Node;
class Myosin_Node;
class Filament;
//=====================
// Include Declarations
#include <iostream>
#include <vector>
#include <memory>
#include "params.h"
#include "coord.h"
#include "externs.h"
//=====================

//==============================
//Parent class:
//==============================
class Node{
    protected:
    //class variables shared by all node types defined in derived classes
        Coord my_location; //every node has a position
        Coord new_total_force; //every node has a total force (sum of spring, bending, stochastic, etc forces) acting on it
        double drag_coeff; //although they have different values, there is an actin damping coeff and a myosin damping coeff
        int my_node_rank; //which node am I; Ex: 1st, 2nd, 3rd in the filament
        
    public:
    //functions that will/can be performed on all types of nodes
    //Constructor:
    //=======================
    Node(Coord location);

    //Destructor
    //=======================
    virtual ~Node();

    //Getters & setters
    //=======================
    virtual Coord get_Node_Location(){return my_location;}
    virtual Coord get_Total_Force(){return new_total_force;}
    virtual double get_Drag_Coeff(){return drag_coeff;}
    virtual void set_Drag_Coeff(double drag);
    virtual int get_My_Node_Rank(){return my_node_rank;}
    virtual void set_My_Node_Rank(int node_num);

    //Functions
    //=======================




};

//==============================
//Derived classes:
//==============================

//Actin class:
//==============================
class Actin_Node: public Node, public enable_shared_from_this<Actin_Node>{
    protected:
    //variables shared by actin nodes
        shared_ptr<Filament> my_filament;
        shared_ptr<Actin_Node> left_neighbor;
        shared_ptr<Actin_Node> right_neighbor;

        //needed for the linear spring force calculation:
        double k_linear_actin;
        double actin_spring_equi_len;
                
        //needed for the bending spring force calculation:
        double k_bend_actin;
        double my_current_angle;
        double equi_angle;

    public:
    //Constructors:
    //=======================
    Actin_Node(Coord location, shared_ptr<Filament> my_filament);
    Actin_Node(Coord location, shared_ptr<Filament> my_filament, shared_ptr<Actin_Node> left_nbh, shared_ptr<Actin_Node> right_nbh);

    //Destructor
    //=======================
    ~Actin_Node();

    //Getters & setters
    //=======================
    shared_ptr<Filament> get_My_Filament(){return my_filament;}
    void set_My_Filament(shared_ptr<Filament> filament);

    shared_ptr<Actin_Node> get_Left_Neighbor(){return left_neighbor;}
    void set_Left_Neighbor(shared_ptr<Actin_Node> left_nbh);

    shared_ptr<Actin_Node> get_Right_Neighbor(){return right_neighbor;}
    void set_Right_Neighbor(shared_ptr<Actin_Node> right_nbh);


    //Functions
    //=======================

};

//Myosin class:
//==============================

//===========================
#endif 