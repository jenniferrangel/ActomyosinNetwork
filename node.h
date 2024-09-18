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
class Myosin;
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
        int vtk_index;    //every node has a vtk ID
        
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
    virtual int get_VTK_Index(){return vtk_index;}
    virtual void update_VTK_Index(int id);

    //Functions
    //=======================
    //Update locations via Langevin Equation
    virtual void update_Position(int Ti);



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

        //for the actomyosin connection
        shared_ptr<Myosin_Node> curr_conn_myo_node;
        shared_ptr<Myosin> curr_conn_minifilament;

        //needed for the linear spring force calculation:
        double k_linear_actin;
        double actin_spring_equi_len;
                
        //needed for the bending spring force calculation:
        double k_bend_actin;
        double my_current_angle;
        double equi_angle;
        double cross_product;  //this is used to determine if angle > 180

        //needed for the actomyosin conn linear spring force calculation:
        double k_linear_actomyoConn;
        double actomyo_spring_equi_len;

        //myo pulling force
        double myo_pull_force;

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

    //set & get the left and right neighbors:
    shared_ptr<Actin_Node> get_Left_Neighbor(){return left_neighbor;}
    void set_Left_Neighbor(shared_ptr<Actin_Node> left_nbh);

    shared_ptr<Actin_Node> get_Right_Neighbor(){return right_neighbor;}
    void set_Right_Neighbor(shared_ptr<Actin_Node> right_nbh);

    //for the actomyosin connection
    shared_ptr<Myosin_Node> get_Conn_Myosin_Node(){return curr_conn_myo_node;}
    void set_Connected_Myosin_Node(shared_ptr<Myosin_Node> myo_node);

    shared_ptr<Myosin> get_Conn_Myosin_MiniFilament(){return curr_conn_minifilament;}
    void set_Conn_Myosin_MiniFilament(shared_ptr<Myosin> minifilament);

    //set & get the linear spring constant and equilibrium length:
    double get_K_Linear_Actin(){return k_linear_actin;}
    void set_K_Linear_Actin(double k_linear);
    double get_Actin_Equi_Len(){return actin_spring_equi_len;}
    void set_Actin_Equi_Len(double equi_len);

    //set & get k_bend, current angle and equilibrium angle:
    double get_K_Bend_Actin(){return k_bend_actin;}
    void set_K_Bend_Actin(double k_bend);

    double get_Current_Angle(){return my_current_angle;}
    void calc_Current_Angle();
    void set_Current_Angle(double my_angle);

    double get_Equi_Angle(){return equi_angle;}
    void set_Equi_Angle(double new_angle);

    //set & get spring constant and equilibrium length for actomyosin connection:
    double get_K_Linear_Actomyo_Conn(){return k_linear_actomyoConn;}
    void set_K_Linear_Actomyo_Conn(double k_lin);
    
    double get_Actomyo_Equi_Len(){return actomyo_spring_equi_len;}
    void set_Actomyo_Equi_Len(double equil_len);

    //set & get myosin pulling force
    double get_Myo_Pull_Force(){return myo_pull_force;}
    void set_Myo_Pull_Force(double pull_force);
    

    //Functions
    //=======================
    void sound_Off_Node_Info();
    void sound_Off_Neighbors();

    //***Functions for calculating forces
    void calculate_Forces(int Ti);

    //Linear spring forces
    Coord calc_Linear_Force();
    Coord linear_Spring_Equation(shared_ptr<Actin_Node> node);

    //Bending spring force at every node triplet
    Coord calc_Bending_Force();
    Coord bending_Force_Equation_Center();
    Coord bending_Force_Equation_Left();
    Coord bending_Force_Equation_Right();

    //Random force
    Coord calc_Stochastic_Force();
    double get_Random_Number(double mean, double stddev);

    //Actomyosin connection force
    Coord calc_Actomyo_Conn_Force();
    Coord linear_Actomyo_Spring_Equation(shared_ptr<Myosin_Node> myo_node);

    //Myosin pulling force from conn
    Coord calc_Myosin_Pull_Force();
    Coord myosin_Pulling_Force(shared_ptr<Myosin_Node> myo_node);


};

//Myosin class:
//==============================
class Myosin_Node: public Node, public enable_shared_from_this<Myosin_Node>{
    protected:
    //variables shared by myosin nodes
    shared_ptr<Myosin> my_myosin_minifilament;
    shared_ptr<Myosin_Node> neighboring_pair;   //every mini-filament is composed of 2 nodes, so each node has it's pair i.e. its neighbor

    //For actomyosin connection
    bool connected;   //false/0: not connected true/1:connected
    shared_ptr<Filament> curr_conn_filament;     // current actin filament it is connected to
    shared_ptr<Actin_Node> curr_conn_actin_node; // current actin node it binds to to form connection
    shared_ptr<Filament> prev_conn_filament;
    shared_ptr<Actin_Node> prev_conn_actin_node;

    //needed for the linear spring force calculation:
    double k_linear_myosin;
    double myosin_spring_equi_len;

    //needed for the actomyosin conn linear spring force calculation:
    double k_linear_actomyo_conn;
    double actomyo_spring_equi_length;

    //myo pulling force
    double myo_pulling_force;

    public:
    //Constructors:
    //=======================
    Myosin_Node(Coord location, shared_ptr<Myosin> my_myosin_minifilament);
    Myosin_Node(Coord location, shared_ptr<Myosin> my_myosin_minifilament, shared_ptr<Myosin_Node> nbring_pair);

    //Destructor
    //=======================
    ~Myosin_Node();

    //Getters & setters
    //=======================
    shared_ptr<Myosin> get_My_Myosin(){return my_myosin_minifilament;}
    void set_My_Myosin(shared_ptr<Myosin> minifilament);

    //set & get the neighboring myosin node pair
    shared_ptr<Myosin_Node> get_Neighboring_Pair(){return neighboring_pair;}
    void set_Neighboring_Pair(shared_ptr<Myosin_Node> nbh_pair);

    //set & get if myosin is connected
    bool get_If_Connected(){return connected;}
    void set_If_Connected(bool conn);

    shared_ptr<Filament> get_Connected_Filament(){return curr_conn_filament;}
    void set_Connected_Filament(shared_ptr<Filament> filament);

    shared_ptr<Actin_Node> get_Connected_Actin_Node(){return curr_conn_actin_node;}
    void set_Connected_Actin_Node(shared_ptr<Actin_Node> node);

    shared_ptr<Filament> get_Prev_Connected_Filament(){return prev_conn_filament;}
    void set_Prev_Connected_Filament(shared_ptr<Filament> prev_filament);

    shared_ptr<Actin_Node> get_Prev_Connected_Actin_Node(){return prev_conn_actin_node;}
    void set_Prev_Connected_Actin_Node(shared_ptr<Actin_Node> prev_node);

    //set & get the linear spring constant and equilibrium length:
    double get_K_Linear_Myosin(){return k_linear_myosin;}
    void set_K_Linear_Myosin(double k_linear);
    double get_Myosin_Equi_Len(){return myosin_spring_equi_len;}
    void set_Myosin_Equi_Len(double equi_len);

    //set & get the spring constant and equilibrium length for the actomyosin connection:
    double get_K_Linear_Actomyo_Conn(){return k_linear_actomyo_conn;}
    void set_K_Linear_Actomyo_Conn(double k_actomyo);

    double get_Actomyo_Spring_Equi_Length(){return actomyo_spring_equi_length;}
    void set_Actomyo_Spring_Equi_Length(double equi_len);

    //set & get myo pulling force
    double get_Myo_Pulling_Force(){return myo_pulling_force;}
    void set_Myo_Pulling_Force(double pulling_force);

    //Functions
    //=======================
    void sound_Off_Nyosin_Node_Info();
    void sound_Off_Neighboring_Pair();

    //***Functions for calculating forces
    void calculate_Myosin_Forces(int Ti);

    //Linear spring forces
    Coord calc_Myosin_Linear_Force();
    Coord linear_Myosin_Spring_Equation(shared_ptr<Myosin_Node> node);

    //Random force
    Coord calc_Myosin_Stochastic_Force();
    double get_Random_Num(double mean, double stddev);

    //Form actomyosin connections:
    void connect_To_Filament(double radius);
    Coord calc_Actomyo_Connection_Force();
    Coord linear_Spring_Actomyo_Equation(shared_ptr<Actin_Node> actin_node);

    //Myosin pulling force from connection
    Coord calc_Myosin_Pulling_Force();
    Coord myosin_Pulling_Force(shared_ptr<Filament> conn_filament, shared_ptr<Actin_Node> actin_node);

};


//===========================
#endif 