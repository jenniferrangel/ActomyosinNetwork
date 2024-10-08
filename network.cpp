//network.cpp
//=========================

//=========================
//Include Dependencies
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <memory>
#include <algorithm>
#include <random>
#include "params.h"
#include "coord.h"
#include "filament.h"
#include "myosin.h"
#include "network.h"
//=========================
// Public Member Functions for network.cpp

//=======================================================
// Constructors
//=======================================================
//network constructor makes new network from .txt file that is read in

/* Network::Network(string filename, mt19937 gen){
    num_filaments = 0;
    this->gen = gen;

    //open the file
    ifstream ifs(filename.c_str());
    
    //check if input file is open
    if(!ifs.is_open()) {
		cout << filename << " is not available! Please check input file name." << endl;
		return;
	}

    stringstream ss;
    string line;
    string temp;
    char garbage;

    //variables to initialize a filament
    int fil_num;
    bool firstIsBarbed;
    bool lastIsBarbed;
    Coord initialNode;
    double x;
    double y;
    Network* my_network = this;

    //parse through the input file line by line and read into "line"
    while(getline(ifs,line)){
        ss.str(line);

        getline(ss,temp,':'); //read values separated by : and stored into "temp"

        if (temp == "FilamentNumber"){
            ss >> fil_num;
        }
        else if(temp == "FirstIsBardedEnd"){
            ss >> firstIsBarbed;            
        }
        else if(temp == "LastIsBarbedEnd"){
            ss >> lastIsBarbed;
        }
        else if(temp == "Node1"){
            ss >> x >> garbage >> y;
            Coord location(x,y);
            initialNode = location;            
        }
        else if(temp == "End_Filament"){
            //input data is used to create a new filament. New filament is pushed onto vector that holds all filaments in the network.
            cout << "Making an actin filament" << endl;

            shared_ptr<Filament> curr = make_shared<Filament>(my_network, fil_num, firstIsBarbed, lastIsBarbed, initialNode);
            //give the filaments its nodes:
            curr->make_nodes(); 
            num_filaments++;
            //push filament onto vector that holds all actin filaments in tissue
            filaments.push_back(curr); 
        }
        ss.clear();
    }
    ifs.close();

    //The following is for checking purposes:
    cout << "My number of filaments is " << get_Num_Filaments() << endl;
    // cout << "My filament # is: " << fil_num << endl;
    // if(firstIsBarbed == true){
    //     cout << "First node is barbed end: " << firstIsBarbed << endl;
    // }
    // cout << "Last node is barbed end: " << lastIsBarbed << endl;
    // cout << "My location is: " << initialNode << endl;
} */

Network::Network(string filename, mt19937 gen, bool predefined_nodes){
    if(predefined_nodes){
        num_filaments = 0;
        this->gen = gen;

        //open the file
        ifstream ifs(filename.c_str());
        
        //check if input file is open
        if(!ifs.is_open()) {
            cout << filename << " is not available! Please check input file name." << endl;
            return;
        }

        stringstream ss;
        string line;
        string temp;
        char garbage;

        //variables to initialize a actin filament and myosin mini-filament
        int fil_num;
        bool firstIsBarbed;
        bool lastIsBarbed;
        vector<Coord> nodes;
        double x;
        double y;
        Network* my_network = this;

        int myo_num;
        vector<Coord> myoNodes;
        double x_myo;
        double y_myo;

        //parse through the input file line by line and read into "line"
        while(getline(ifs,line)){
            ss.str(line);

            getline(ss,temp,':'); //read values separated by : and stored into "temp"

            if (temp == "FilamentNumber"){
                ss >> fil_num;
                cout << "Filament num = " << fil_num << endl;
            }
            else if(temp == "FirstIsBardedEnd"){
                ss >> firstIsBarbed;
                cout << firstIsBarbed << endl;            
            }
            else if(temp == "LastIsBarbedEnd"){
                ss >> lastIsBarbed;
                cout << lastIsBarbed << endl;
            }
            else if(temp == "Node"){
                cout << "Saving node position into vector" << endl;
                ss >> x >> garbage >> y;
                cout << x << "," << y << endl;
                Coord location(x,y);
                nodes.push_back(location);
                cout << location.get_X() << "," << location.get_Y() << endl;
            }
            else if(temp == "End_Filament"){
                //input data is used to create a new filament. New filament is pushed onto vector that holds all filaments in the network.
                cout << "Making an actin filament..." << endl;

                shared_ptr<Filament> curr = make_shared<Filament>(my_network, fil_num, firstIsBarbed, lastIsBarbed);
                //give the filaments its nodes:
                curr->make_nodes(nodes); 
                num_filaments++;
                //push filament onto vector that holds all actin filaments in tissue
                filaments.push_back(curr); 

                //to avoid carrying nodes from other filaments reset the nodes vector
                nodes.clear();
                cout << "Actin filament created!" << endl;
            }
            else if(temp == "MyosinNumber"){
                ss >> myo_num;  
                cout << "MyosinNumber = " << myo_num << endl;          
            }
            else if(temp == "MyoNode"){
                cout << "Saving myosin node position into vector" << endl;
                ss >> x_myo >> garbage >> y_myo;
                cout << x_myo << "," << y_myo << endl;
                Coord loc(x_myo,y_myo);
                myoNodes.push_back(loc);
                cout << loc.get_X() << "," << loc.get_Y() << endl;
            }
            else if(temp == "End_Myosin"){
                //input data is used to create a new myosin mini-filament. New myosin is pushed onto vector that holds all myosin minifilaments in the network.
                cout << "Making a myosin mini-filament..." << endl;

                shared_ptr<Myosin> curr = make_shared<Myosin>(my_network, myo_num);
                //give the myosin mini-filaments its nodes:
                curr->make_myosin_nodes(myoNodes);
                num_myosins++;
                //push myosin mini-filament onto vector that holds all myosins in tissue
                myosins.push_back(curr);

                //to avoid carrying nodes from other myosins reset the nodes vector
                myoNodes.clear();
                cout << "Myosin mini-filament created!" << endl;
            }
            ss.clear();
        }
        ifs.close();

        //The following is for checking purposes:
        cout << "My number of filaments is " << get_Num_Filaments() << endl;
        // cout << "My filament # is: " << fil_num << endl;
        // if(firstIsBarbed == true){
        //     cout << "First node is barbed end: " << firstIsBarbed << endl;
        // }
        // cout << "Last node is barbed end: " << lastIsBarbed << endl;
        cout << "My number of myosin mini-filaments is " << get_Num_Myosins() << endl;

    }else{
        cout << "Need to create the nodes in the network" << endl;

        num_filaments = 0;
        this->gen = gen;

        //open the file
        ifstream ifs(filename.c_str());
        
        //check if input file is open
        if(!ifs.is_open()) {
            cout << filename << " is not available! Please check input file name." << endl;
            return;
        }

        stringstream ss;
        string line;
        string temp;
        char garbage;

        //variables to initialize a filament
        int fil_num;
        bool firstIsBarbed;
        bool lastIsBarbed;
        Coord initialNode;
        double x;
        double y;
        Network* my_network = this;

        //parse through the input file line by line and read into "line"
        while(getline(ifs,line)){
            ss.str(line);

            getline(ss,temp,':'); //read values separated by : and stored into "temp"

            if (temp == "FilamentNumber"){
                ss >> fil_num;
            }
            else if(temp == "FirstIsBardedEnd"){
                ss >> firstIsBarbed;            
            }
            else if(temp == "LastIsBarbedEnd"){
                ss >> lastIsBarbed;
            }
            else if(temp == "Node1"){
                ss >> x >> garbage >> y;
                Coord location(x,y);
                initialNode = location;            
            }
            else if(temp == "End_Filament"){
                //input data is used to create a new filament. New filament is pushed onto vector that holds all filaments in the network.
                cout << "Making an actin filament..." << endl;

                shared_ptr<Filament> curr = make_shared<Filament>(my_network, fil_num, firstIsBarbed, lastIsBarbed, initialNode);
                //give the filaments its nodes:
                curr->make_nodes(); 
                num_filaments++;
                //push filament onto vector that holds all actin filaments in tissue
                filaments.push_back(curr); 

                cout << "Actin filament created!" << endl;
            }
            ss.clear();
        }
        ifs.close();

        //The following is for checking purposes:
        cout << "My number of filaments is " << get_Num_Filaments() << endl;
        // cout << "My filament # is: " << fil_num << endl;
        // if(firstIsBarbed == true){
        //     cout << "First node is barbed end: " << firstIsBarbed << endl;
        // }
        // cout << "Last node is barbed end: " << lastIsBarbed << endl;
        // cout << "My location is: " << initialNode << endl;
        }
}

//=======================================================
// Getters and Setters
//=======================================================
//returning or updating the number of actin filaments in the network:
void Network::get_Filaments(vector<shared_ptr<Filament>>& filaments){
    filaments = this->filaments;
    return;
}

void Network::update_Num_Filaments(shared_ptr<Filament>& new_Filament){
   num_filaments++;
   filaments.push_back(new_Filament); 
   return;
}

//Generate and return a random number from a normal distribution based on the specified mean and standard deviation
double Network::get_Normally_Distributed_Random_Number(double mean, double stddev){
    std::normal_distribution<double> distribution(mean,stddev);
	double random_num = distribution(this->gen);
	return random_num;
}

void Network::get_Myosins(vector<shared_ptr<Myosin>>& myosins){
    myosins = this->myosins;
    return;
}

void Network::update_Num_Myosins(shared_ptr<Myosin>& new_Myosin){
    num_myosins++;
    myosins.push_back(new_Myosin);
    return;
}

//=======================================================
// Functions
//=======================================================
void Network::sound_Off_All_Node_Info(){
    cout << "********************************************************************************************************************************" << endl;
    cout << "The network is composed of " << get_Num_Filaments() << " actin filaments AND " << get_Num_Myosins() << " myosin mini-filaments!!" << endl;
    cout << "********************************************************************************************************************************" << endl;

    cout << "Actin filament info: " << endl;
    cout << "====================" << endl;

    shared_ptr<Filament> curr;
    vector<shared_ptr<Actin_Node>> actins;

    for(unsigned int i = 0; i < filaments.size(); i++){
        curr = filaments.at(i);
        cout << "Filament #: " << curr->get_Filament_Num()<< endl;
        cout << "First node is barbed: " << curr->get_First_Node_Polarity()<< endl;
        cout << "Last node is barbed: " << curr->get_Last_Node_Polarity()<< endl;
        cout << "Total number of nodes in filament: " << curr->get_Num_Actin_Nodes()<< endl;

        curr->get_Actin_Nodes_Vec(actins);

        for(unsigned int j = 0; j < actins.size(); j++){
            actins.at(j)->sound_Off_Node_Info();
        }
        cout << endl;
        cout << "***********************" << endl;
        cout << endl;
    }

    cout << "Myosin mini-filament info: " << endl;
    cout << "==========================" << endl;

    shared_ptr<Myosin> current;
    vector<shared_ptr<Myosin_Node>> minifilaments;
    for(unsigned int i = 0; i < myosins.size(); i++){
        current = myosins.at(i);
        cout << "Myosin mini-filament #: " << current->get_Myosin_Num() << endl;
        cout << "Total number of nodes in mini-filament: " << current->get_Num_Myosin_Nodes() << endl;

        current->get_Myosin_Nodes_Vec(minifilaments);

        for(unsigned int j = 0; j < minifilaments.size(); j++){
            minifilaments.at(j)->sound_Off_Nyosin_Node_Info();
        }
        cout << endl;
        cout << "***********************" << endl;
        cout << endl;
    }

    cout << "********************************************************************************************************************************" << endl;
    
    return;
}

void Network::sound_Off_Neighbors(){
    cout << "ACTIN NEIGHBORS: " << endl;
    cout << "***********************************************" << endl;

    shared_ptr<Filament> curr;
    vector<shared_ptr<Actin_Node>> actins;

    for(unsigned int i = 0; i < filaments.size(); i++){
        curr = filaments.at(i);
        cout << "Filament #: " << curr->get_Filament_Num()<< endl;

        curr->get_Actin_Nodes_Vec(actins);

        for(unsigned int j = 0; j < actins.size(); j++){
            actins.at(j)->sound_Off_Neighbors();
        }
    }

    cout << "***********************************************" << endl;

    cout << endl;

    cout << "MYOSIN NEIGHBORING PAIRS: " << endl;
    cout << "***********************************************" << endl;

    shared_ptr<Myosin> current;
    vector<shared_ptr<Myosin_Node>> minifilaments;

    for(unsigned int i = 0; i < myosins.size(); i++){
        current = myosins.at(i);
        cout << "Myosin mini-filament #: " << current->get_Myosin_Num() << endl;
        
        current->get_Myosin_Nodes_Vec(minifilaments);

        for(unsigned int j = 0; j < minifilaments.size(); j++){
            minifilaments.at(j)->sound_Off_Neighboring_Pair();
        }
    }

    cout << "***********************************************" << endl;

    return;
}

void Network::sound_Off_Possible_Connections(int Ti){
    cout << "POSSIBLE CONNECTIONS: " << endl;
    cout << "***********************************************" << endl;

    cout << "Timestep = " << Ti << endl;

    shared_ptr<Myosin> current;
    // vector<shared_ptr<Myosin_Node>> myosinNodes;
    vector<shared_ptr<Filament>> actinfilaments;

    for(unsigned int i = 0; i < myosins.size(); i++){
        current = myosins.at(i);
        cout << "Myosin mini-filament #: " << current->get_Myosin_Num() << endl;
        cout << "My neighboring actin filaments are: " << endl;
        actinfilaments = current -> get_Possible_Connections();
        if(actinfilaments.empty()){
            cout << "NO CONNECTIONS" << endl;
        }
        else{
            for(unsigned int j = 0; j < actinfilaments.size(); j++){
                cout << actinfilaments.at(j)->get_Filament_Num() << endl;
            }
        }
    }
}

//***Calculates forces acting on each filament node***//
void Network::calculate_New_Forces(int Ti){
   //#pragma omp parallel for schedule(static,1)
   for(unsigned int i = 0; i < filaments.size(); i++){
        cout << "***********************************************" << endl;
        cout << "ACTIN FILAMENT FORCES BEING CALCULATED..." << endl;
        cout << "***********************************************" << endl;
        cout << "Calc forces for actin filament: " << i << endl;
        filaments.at(i)->calculate_New_Forces(Ti);
        cout << "Successfully calculated actin forces!!" << endl;
    }

    //#pragma omp parallel for schedule(static,1)
   for(unsigned int j = 0; j < myosins.size(); j++){
        cout << "***********************************************" << endl;
        cout << "MYOSIN MINI-FILAMENT FORCES BEING CALCULATED..." << endl;
        cout << "***********************************************" << endl;
        cout << "Calc forces for myosin mini-filament: " << j << endl;
        myosins.at(j)->calculate_New_Myosin_Forces(Ti);
        cout << "Successfully calculated myosin forces!!" << endl;
    }
   
   return;

}

//***Update node positions via Langevin equation***//
void Network::update_Positions(int Ti){
    //#pragma omp parallel for schedule(static,1)
    for(unsigned int i = 0; i < filaments.size(); i++){
        cout << "UPDATING ACTIN NODE POSITIONS..." << endl;
        cout << "***********************************************" << endl;
        filaments.at(i)->update_Node_Positions(Ti);
    }

    //#pragma omp parallel for schedule(static,1)
    for(unsigned int j = 0; j < myosins.size(); j++){
        cout << "UPDATING MYOSIN NODE POSITIONS..." << endl;
        cout << "***********************************************" << endl;
        myosins.at(j)->update_Myosin_Node_Positions(Ti);
    }

    return;
}

//***Myosin searches through each filament to find filaments within reach that it can form connections with***//
void Network::find_Possible_Connections(int Ti){
    double conn_radius = ACTIN_MYO_CONNECT_RADIUS;
    cout << "connection radius = " << conn_radius << endl;

    //#pragma omp parallel for schedule(static,1)
    for(unsigned int i = 0; i < myosins.size(); i++){
        cout << "MYOSIN SEARCHING FOR NEIGHBORING ACTIN FILAMENTS..." << endl;
        cout << "***********************************************" << endl;
        myosins.at(i)->find_Possible_FilConn(filaments, conn_radius);
    }
    return;
}

void Network::update_Possible_Connections(int Ti){
    double conn_radius = ACTIN_MYO_CONNECT_RADIUS;
    cout << "connection radius = " << conn_radius << endl;

    //empty out old connections
    for(unsigned int i = 0; i < myosins.size(); i++){
        myosins.at(i)->clear_Possible_Connections();
    }

    //find new ones
    //#pragma omp parallel for schedule(static,1)
    for(unsigned int j = 0; j < myosins.size(); j++){
        cout << "UPDATING..MYOSIN SEARCHING FOR NEIGHBORING ACTIN FILAMENTS..." << endl;
        cout << "***********************************************" << endl;
        myosins.at(j)->find_Possible_FilConn(filaments, conn_radius);
    }

    return;
}

//***Form actomyosin connections***//
void Network::formActomyoConnections(int Ti){
    //Important notes:
    //1)Myosin will ALWAYS move towards barbed end
    //2)The 2 myosin heads CANNOT bind to the same actin

    double conn_radius = ACTIN_MYO_CONNECT_RADIUS;
    cout << "connection radius = " << conn_radius << endl;

    for(unsigned int i = 0; i < myosins.size(); i++){
        cout << "FORMING ACTOMYOSIN CONNECTIONS..." << endl;
        cout << "***********************************************" << endl;
        myosins.at(i)->formActomyoConnections(Ti, conn_radius);
    }
    
    return;
}

//***Functions for VTK output***//
//Actin Filaments
void Network::print_VTK_File(ofstream& ofs){

    //Update the index of the nodes
    update_VTK_Indices();

    ofs << "# vtk DataFile Version 3.0" << endl;
	ofs << "Point representing Actin Filaments in the Actomyosin Network model" << endl;
	ofs << "ASCII" << endl << endl;
	ofs << "DATASET UNSTRUCTURED_GRID" << endl;

    //Set the structure:
    //Need the total number of points/nodes for all filaments and number of edges/connections per filament
    int num_Points = 0;
    int num_Edges = 0; 

    for(unsigned int i = 0; i < filaments.size(); i++){
        num_Points += filaments.at(i)->get_Num_Actin_Nodes();
        num_Edges += filaments.at(i)->get_Num_Actin_Nodes() - 1; //always 1 less than total # nodes per filament
        cout << "number of edges = " << num_Edges << endl;
    }

    //The node/point positions
    ofs << "POINTS " << num_Points << " float64" << endl;

    // vector<int> start_points;
	// vector<int> end_points;
    //for double-checking purposes so we know it is getting a total # nodes
	int count = 0;

    for(unsigned int i = 0; i < filaments.size(); i++){
        //start_points.push_back(count);
        filaments.at(i)->print_VTK_Points(ofs,count);
        //end_points.push_back(count - 1);
    }
    //cout << "count = " << count << endl;

    ofs << endl;

    //to visualize how the nodes are connected i.e. the springs connecting nodes:
    ofs << "CELLS " << num_Edges << " " << 3*num_Edges << endl;
    for(unsigned int i = 0; i < filaments.size(); i++){
        filaments.at(i)->print_VTK_connections(ofs);
    }

    ofs << endl;

    ofs << "CELL_TYPES " << num_Edges << endl;
    for(int i = 0; i < num_Edges; i++){
        //type for adhesion relationship: in this case a vtk_line ._____.
        ofs << 3 << endl;
    }

    ofs << endl;

    //Set the attributes:
    ofs << "POINT_DATA " << num_Points << endl;
    //To visualized the barbed end
    ofs << "SCALARS Barbed_End int " << 1 << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for(unsigned int i = 0; i < filaments.size(); i++){
        filaments.at(i)->print_VTK_BarbedEnds(ofs);
    }

    ofs << endl;

    return;

}

void Network::update_VTK_Indices(){
    cout << "Updating Actin VTK indices..." << endl;

    int id = 0;
    //iterate through the filaments to reassign vtk id's (starting at 0)
    for (unsigned int i = 0; i < filaments.size(); i++){
        filaments.at(i)->update_Node_VTK_Indices(id);
        //cout << "ID = " << id << endl;
    }

    //cout << "final index = " << id << endl;

    return;

}

//Myosin Mini-filaments
void Network::print_Myosin_VTK_File(ofstream& ofs){
    //Update the index of the nodes
    update_Myosin_VTK_Indices();

    ofs << "# vtk DataFile Version 3.0" << endl;
	ofs << "Point representing Myosin MiniFilaments in the Actomyosin Network model" << endl;
	ofs << "ASCII" << endl << endl;
	ofs << "DATASET UNSTRUCTURED_GRID" << endl;

    //Set the structure:
    //Need the total number of points/nodes for all myosin mini-filaments and number of edges/connections per mini-filament
    int num_Points_Myo = 0;
    int num_Edges_Myo = 0;

    for(unsigned int i = 0; i < myosins.size(); i++){
        num_Points_Myo += myosins.at(i)->get_Num_Myosin_Nodes();
        num_Edges_Myo += myosins.at(i)->get_Num_Myosin_Nodes() - 1; //always 1 less than total # nodes per mini-filament
        cout << "number of myosin edges = " << num_Edges_Myo << endl;
    }

    //The myosin node/point positions
    ofs << "POINTS " << num_Points_Myo << " float64" << endl;

    // vector<int> start_points;
	// vector<int> end_points;
    //for double-checking purposes so we know it is getting a total # nodes
	int count_myo = 0;

    for(unsigned int i = 0; i < myosins.size(); i++){
        //start_points.push_back(count_myo);
        myosins.at(i)->print_Myosin_VTK_Points(ofs,count_myo);
        //end_points.push_back(count_myo - 1);
    }
    cout << "count = " << count_myo << endl;

    ofs << endl;

    //to visualize how the nodes are connected i.e. the springs connecting nodes:
    ofs << "CELLS " << num_Edges_Myo << " " << 3*num_Edges_Myo << endl;
    for(unsigned int i = 0; i < myosins.size(); i++){
        myosins.at(i)->print_Myosin_VTK_connections(ofs);
    }

    ofs << endl;

    ofs << "CELL_TYPES " << num_Edges_Myo << endl;
    for(int i = 0; i < num_Edges_Myo; i++){
        //type for adhesion relationship: in this case a vtk_line ._____.
        ofs << 3 << endl;
    }

    ofs << endl;


    return;
}

void Network::update_Myosin_VTK_Indices(){
    cout << "Updating Myosin VTK indices..." << endl;

    int id = 0;
    //iterate through the myosin mini-filaments to reassign vtk id's (starting at 0)
    for (unsigned int i = 0; i < myosins.size(); i++){
        myosins.at(i)-> update_Myosin_Node_VTK_Indices(id);
        //cout << "ID = " << id << endl;
    }
    //cout << "final index = " << id << endl;

    return;
}

//***Functions for printing data output***//
//This function prints out the node locations
void Network::locations_Output(ofstream& ofs, int Ti){
    ofs <<"Filament#" << ' ' << "Node#" << ' ' << "Location (x, y, z)" << endl;
    for(unsigned int i = 0; i < filaments.size(); i++){
        filaments.at(i)->print_Location(ofs, Ti);
    }

    return;
}

//This function prints the actin node data info
void Network::node_Data_Output(ofstream& ofs, int Ti){
    ofs <<"Filament#" << ' ' << "Node#" << ' ' << "Location (x, y, z)" << ' ' <<  "Left Nbr (x, y, z)" << ' ' <<  "Right Nbr (x, y, z)" << ' ' << "Drag Coeff" << ' ' << "k_linear" << ' ' << "Equi Len" << ' ' << "k_bend" << ' ' << "theta" << ' ' << "theta_equi" << endl;
    for(unsigned int i = 0; i < filaments.size(); i++){
        filaments.at(i)->print_Node_Data(ofs, Ti);
    }
    return;
}

//This function prints the filament data node info
void Network::filament_Data_Output(ofstream& ofs, int Ti){
    ofs << "Filament#" << ' ' << "#Nodes" << ' ' << "Polarity_1stNode" << ' ' << "Polarity_LastNode" << endl;
    for(unsigned int i = 0; i < filaments.size(); i++){
        filaments.at(i)->print_Filament_Data(ofs, Ti);
    }
    return;
}

//This function prints the network data node info
void Network::network_Data_Output(ofstream& ofs, int Ti){
    ofs << "NumberActinFils" << ' ' << "NumberMyosins" << ' ' << "Ti" << endl;
    ofs << this->get_Num_Filaments() << ' ' << this->get_Num_Myosins() << ' ' << Ti << endl;
    return;
}

//This function prints the myosin node data info
void Network::locations_Myosin_Output(ofstream& ofs, int Ti){
    ofs <<"Myosin#" << ' ' << "Node#" << ' ' << "Location (x, y, z)" << endl;
    for(unsigned int i = 0; i < myosins.size(); i++){
        myosins.at(i)->print_Myosin_Locations(ofs, Ti);
    }
    
    return;
}

//This function prints the myosin node data info
void Network::myosin_Node_Data_Output(ofstream& ofs, int Ti){
    ofs <<"Myosin#" << ' ' << "Node#" << ' ' << "Location (x, y, z)" << ' ' <<  "Nbr Pair (x, y, z)" << ' ' << "Drag Coeff" << ' ' << "k_linear" << ' ' << "Equi Len" << endl;
    for(unsigned int i = 0; i < myosins.size(); i++){
        myosins.at(i)->print_Myosin_Node_Data(ofs, Ti);
    }
    return;
}

//This function prints the myosin mini-filament data info
void Network::Myosin_Minifilament_Data_Output(ofstream& ofs, int Ti){
    ofs << "Myosin#" << ' ' << "#Nodes" << endl;
    for(unsigned int i = 0; i < myosins.size(); i++){
        myosins.at(i)->print_MiniFilament_Data(ofs, Ti);
    }
    return;
}

//=======================================================
// Destructor
//=======================================================
Network::~Network() {
	cout << "Network destructor executed!" << endl;
}
//=========================
//End of network.cpp

