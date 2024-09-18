//myosin.cpp:
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
#include "myosin.h"
#include "network.h"
#include <boost/random.hpp>
//===================

//==========================================================
/** Myosin Mini-filament Class Member Functions**/

//==============================
// Constructors
//==============================
Myosin::Myosin(Network* network, int myosin_Number){
    cout << "Calling the myosin constructor!" << endl;

    this->my_network = network;
    this->myosinNumber = myosin_Number;

    //set the damping/drag coefficient
    this->myosin_drag_coeff = MYO_DRAG_COEFF;
    cout << "My drag coeff is: " << myosin_drag_coeff << endl;

    num_myosin_nodes = 0;

}

//make nodes with predefined positions
void Myosin::make_myosin_nodes(vector<Coord> nodes){
    cout << "Making the myosin nodes!" << endl;

    double current_x;
    double current_y;
    Coord location;
    Coord node;
    shared_ptr<Myosin> this_myosin = shared_from_this();
    //shared_ptr<Myosin_Node> previousMyosin;
    shared_ptr<Myosin_Node> currentMyosin;

    // shared_ptr<Myosin_Node> first; //to save the starter node
    // shared_ptr<Myosin_Node> second;
    //to fix the neighbors of first and second node
    shared_ptr<Myosin_Node> firstNode;
    shared_ptr<Myosin_Node> secondNode;

    //make the successive actin nodes
    for(unsigned int i = 0; i < nodes.size(); i++){
        node = nodes.at(i);
        current_x = node.get_X();
        current_y = node.get_Y();
        location = Coord(current_x, current_y);
        shared_ptr<Myosin_Node> new_node = make_shared<Myosin_Node>(location, this_myosin);
        currentMyosin = new_node;
        myosin_nodes.push_back(currentMyosin);
        num_myosin_nodes++;

        if(i==0){
            //it is the first node in the pair
            firstNode = currentMyosin;
        }
        else if(i==1){
            //it is the second node in the pair
            secondNode = currentMyosin;
        }
    }

    //Assing the neighboring pair: first's neighboring pair is second and vice versa
    firstNode->set_Neighboring_Pair(secondNode);
    secondNode->set_Neighboring_Pair(firstNode);

    //Now we set the private node member variables
    vector<shared_ptr<Myosin_Node>> myosins; //will store myosin nodes here
    this->get_Myosin_Nodes_Vec(myosins);

    //the drag coefficient is inherited from the myosin mini-filament
    double new_drag_coeff = this->get_Myosin_Drag_Coeff();
    // double spring_equi_len = 0.0;
    // double k_linear = 0.0;
    int rank = 0;
    bool conn = false; //all nodes are initialy not connected

    for(unsigned int i = 0; i < myosins.size(); i++){
        rank = i;
        myosins.at(i)->set_Drag_Coeff(new_drag_coeff);
        myosins.at(i)->set_K_Linear_Myosin(K_LINEAR_STIFF_MYOSIN);
        myosins.at(i)->set_Myosin_Equi_Len(MYOSIN_SPRING_EQUI_LEN);
        myosins.at(i)->set_My_Node_Rank(rank);

        //Actomyosin connections are not present when nodes are initialized
        myosins.at(i)->set_If_Connected(conn);
        myosins.at(i)->set_Connected_Filament(nullptr);
        myosins.at(i)->set_Connected_Actin_Node(nullptr);
        myosins.at(i)->set_Prev_Connected_Filament(nullptr);
        myosins.at(i)->set_Prev_Connected_Actin_Node(nullptr);

        //setting up parameters for connections
        myosins.at(i)->set_K_Linear_Actomyo_Conn(K_LINEAR_STIFF_ACTIN_MYO_CONN);
        myosins.at(i)->set_Actomyo_Spring_Equi_Length(ACTIN_MYO_CONN_SPRING_EQUI_LEN);

        //myosin pulling force
        myosins.at(i)->set_Myo_Pulling_Force(F_MYO_PULLING);

    }

    cout << "Finished making the myosin nodes!" << endl;

    return;
}


//==============================
//Getters and setters:
//==============================
void Myosin::set_Myosin_Num(const int rank){
    this->myosinNumber = rank;
    return;
}

void Myosin::set_Myosin_Drag_Coeff(double new_drag){
    this->myosin_drag_coeff = new_drag;
    return;
}

void Myosin::set_Num_Myosin_Nodes(int number_myosin_nodes){
    this->num_myosin_nodes = number_myosin_nodes;
    return;
}

void Myosin::get_Myosin_Nodes_Vec(vector<shared_ptr<Myosin_Node>>& myosins){
    myosins = this->myosin_nodes;
    return;
}

//empty out possible connections
void Myosin::clear_Possible_Connections(){
    cout << "CLEARING possible connections vector" << endl;
    possible_filament_connections.clear();
    return;
}
//==============================
//Functions:
//==============================

//***Calculate forces***//
void Myosin::calculate_New_Myosin_Forces(int Ti){
    vector<shared_ptr<Myosin_Node>> myosins; //will store myosin nodes here
    this->get_Myosin_Nodes_Vec(myosins);

    cout << "Calculating myosin forces..." << endl;

    //#pragma omp parallel
       // {
    //#pragma omp for schedule(static,1)
            for(unsigned int i = 0; i < myosins.size(); i++){
                myosins.at(i)->calculate_Myosin_Forces(Ti);
            }
        //}

    return;
}

//***Update node positions via Langevin Equation***//
void Myosin::update_Myosin_Node_Positions(int Ti){
    vector<shared_ptr<Myosin_Node>> myominifils; //will store myosin nodes here
    this->get_Myosin_Nodes_Vec(myominifils);

    //#pragma omp parallel 
	    //{	
    //#pragma omp for schedule(static,1)
            for(unsigned int i = 0; i < myominifils.size(); i++){
                cout << "Updating myosin locations..." << endl;
                myominifils.at(i)->update_Position(Ti);
            }
        //}
    
    return;
}

//***Find neighboring actin filaments for possible connections***//
void Myosin::find_Possible_FilConn(const vector<shared_ptr<Filament>>& actinfils, double connection_radius){
    cout << "range radius = " << connection_radius << endl;
    //bool foundConnection = false;
    int vectorSize = actinfils.size();
    cout << "vector size = " << vectorSize << endl;

    //Every myosin node must find its own possible connections
    for(unsigned int i = 0; i < myosin_nodes.size(); i++){
        //cout << "ENTERING LOOP i" << endl;

        Coord myo_location = myosin_nodes.at(i)->get_Node_Location();
        cout << "myosin node = (" << myo_location.get_X() << " , " << myo_location.get_Y() << ")" << endl;

        //Loop through every actin filament
        for(unsigned int j = 0; j < actinfils.size(); j++){
            //cout << "ENTERING LOOP j" << endl;

            int fil_num = actinfils.at(j)->get_Filament_Num();
            cout << "Actin filament #:" << fil_num << endl;

            vector<shared_ptr<Actin_Node>> actins; //will store actin nodes here
            actinfils.at(j)->get_Actin_Nodes_Vec(actins);
            shared_ptr<Actin_Node> actin_fil_node = NULL;

            //Go through every node in the actin filament
            for(unsigned int k = 0; k < actins.size(); k++){
                //cout << "ENTERING LOOP k" << endl;

                actin_fil_node = actins.at(k);
                Coord actin_location = actin_fil_node->get_Node_Location();
                cout << "actin node: (" << actin_location.get_X() << " , " << actin_location.get_Y() << ")" << endl;
                
                //Calculate the distance between the actin nodes and myosin node:
                Coord diff_vec = actin_location - myo_location;
                cout << "diff vec = " << actin_location << " - " << myo_location << " = " << diff_vec << endl;

                double distance = diff_vec.length();
                cout << "distance = " << distance << endl;

                //Check is the distance is less than or equal to connection_radius. If yes, it is a possible actin connection
                if(distance <= connection_radius){
                    cout << "Possible connection found" << endl;
                    possible_filament_connections.push_back(actinfils.at(j));
                    break; //move on to the next actin filament
                }

            }
        }

    }



    return;
}

//***Form actomyosin connections***//
void Myosin::formActomyoConnections(int Ti, double connection_radius){
    cout << "MYO MINI-FIL #" << this->get_Myosin_Num() << endl;

    //Get the neighboring actin filaments
    vector<shared_ptr<Filament>> nbhrActinFilaments;
    nbhrActinFilaments = this -> get_Possible_Connections();

    if(nbhrActinFilaments.empty()){
        cout << "NO neighboring actin filaments. NO actomyo connection formed!!" << endl;
    }
    else{
        //for double checking purposes:
        cout << "Neighboring actin filaments:" << endl;
        for(unsigned int i = 0; i < nbhrActinFilaments.size(); i++){
            cout << nbhrActinFilaments.at(i)->get_Filament_Num() << endl;
        }

        //Every myosin node find/form your connection
        for(unsigned int i = 0; i < myosin_nodes.size(); i++){
            
            Coord myo_location = myosin_nodes.at(i)->get_Node_Location();
            cout << "myosin node = " << myo_location << endl;
            
            myosin_nodes.at(i)->connect_To_Filament(connection_radius);
            
        }


    }
    return;
}

//***Functions for VTK output****//
void Myosin::update_Myosin_Node_VTK_Indices(int& id){
    for(unsigned int i = 0; i < myosin_nodes.size(); i++){
        myosin_nodes.at(i)->update_VTK_Index(id);
        cout << "id = " << id << endl;
        id++;
    }
    return;
}

void Myosin::print_Myosin_VTK_Points(ofstream& ofs, int& count_myo){
    for(unsigned int i = 0; i < myosin_nodes.size(); i++){
        Coord location = myosin_nodes.at(i)->get_Node_Location();
        ofs << location.get_X() << ' ' << location.get_Y() << ' ' << 0 << endl;
        count_myo++;
    }

    return;
}

void Myosin::print_Myosin_VTK_connections(ofstream& ofs){
    int my_ID;
    int pair_ID;

    vector<shared_ptr<Myosin_Node>> myosins; //will store myosin nodes here
    this->get_Myosin_Nodes_Vec(myosins);
    shared_ptr<Myosin_Node> curr_node = NULL;
    shared_ptr<Myosin_Node> nhbr_pair = NULL;
    
    for(unsigned int i = 0; i < myosin_nodes.size(); i++){
        curr_node = myosins.at(i);
        nhbr_pair = myosins.at(i)->get_Neighboring_Pair();
        
        //to avoid printing doubles of the same connection between a pair (Ex: 0 1, 1,0 when it should just be 0 1)
        if((i%2)==0){
            cout << "curr_node = " << curr_node->get_Node_Location() << endl;
            cout << "nhbr_pair = "  << nhbr_pair->get_Node_Location() << endl;

            my_ID = curr_node->get_VTK_Index();
            pair_ID = nhbr_pair->get_VTK_Index();
            ofs.flush();
            ofs << 2 << " " << my_ID << " " << pair_ID << endl;

            cout << "my_myosin_ID = " << my_ID << endl;
            cout << "pair_ID = " << pair_ID << endl;
        }
    }

    return;
}

//***Functions for data output***//
//This function prints out the myosin node locations
void Myosin::print_Myosin_Locations(ofstream& ofs, int Ti){
    for(unsigned int i = 0; i < myosin_nodes.size(); i++){
        Coord location = myosin_nodes.at(i)->get_Node_Location();
        ofs << this->get_Myosin_Num() << ' ' << myosin_nodes.at(i)->get_My_Node_Rank() << ' '<< location.get_X() << ' ' << location.get_Y() << ' ' << 0 << endl;
    }
    return;
}

//This function prints out the node data
void Myosin::print_Myosin_Node_Data(ofstream& ofs, int Ti){
    for(unsigned int i = 0; i < myosin_nodes.size(); i++){
        Coord location = myosin_nodes.at(i)->get_Node_Location();
        Coord nbr_pair = myosin_nodes.at(i)->get_Neighboring_Pair()->get_Node_Location();
        
        ofs << this->get_Myosin_Num() << ' ' << myosin_nodes.at(i)->get_My_Node_Rank() << ' '<< location.get_X() << ' ' << location.get_Y() << ' ' << 0 << ' ' << nbr_pair.get_X() << ' ' << nbr_pair.get_Y() << ' ' << 0 << ' ' << myosin_nodes.at(i)->get_Drag_Coeff() << ' ' << myosin_nodes.at(i)->get_K_Linear_Myosin() << ' ' << myosin_nodes.at(i)->get_Myosin_Equi_Len() << endl;
    }

    return;
}

//This function prints out the myosin mini-filament data
void Myosin::print_MiniFilament_Data(ofstream& ofs, int Ti){
    ofs << this->get_Myosin_Num() << ' ' << this->get_Num_Myosin_Nodes() << endl;

    vector<shared_ptr<Filament>> actinfilaments;
    actinfilaments = this -> get_Possible_Connections();
        if(actinfilaments.empty()){
            ofs << "NO POSSIBLE CONNECTIONS" << endl;
        }
        else{
            for(unsigned int j = 0; j < actinfilaments.size(); j++){
                ofs << "Possible actin connections:" << endl;
                ofs << actinfilaments.at(j)->get_Filament_Num() << endl;
            }
        }
    return;
}


//==============================
//Destructor:
//==============================
Myosin::~Myosin(){
    cout << "Myosin destructor executed!" << endl;
}

//==========================================================
// End of myosin.cpp