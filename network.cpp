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
#include "network.h"
//=========================
// Public Member Functions for network.cpp

//=======================================================
// Constructors
//=======================================================
//network constructor makes new network from .txt file that is read in

Network::Network(string filename, mt19937 gen){
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
    cout << "My filament # is: " << fil_num << endl;
    if(firstIsBarbed == true){
        cout << "First node is barbed end: " << firstIsBarbed << endl;
    }
    cout << "Last node is barbed end: " << lastIsBarbed << endl;
    cout << "My location is: " << initialNode << endl;
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

//=======================================================
// Functions
//=======================================================
void Network::sound_Off_All_Node_Info(){
    cout << "The network is composed of " << get_Num_Filaments() << " actin filaments, each of which has " << INIT_NUM_ACTIN_NODES << " nodes." << endl;
    cout << "****************************************************" << endl;

    shared_ptr<Filament> curr;
    vector<shared_ptr<Actin_Node>> actins;

    for(unsigned int i = 0; i < filaments.size(); i++){
        curr = filaments.at(i);
        cout << "Filament #: " << curr->get_Filament_Num()<< endl;
        cout << "First node is barbed: " << curr->get_First_Node_Polarity()<< endl;
        cout << "Last node is barbed: " << curr->get_Last_Node_Polarity()<< endl;
        cout << "Total number of nodes: " << curr->get_Num_Actin_Nodes()<< endl;

        curr->get_Actin_Nodes_Vec(actins);

        for(unsigned int j = 0; j < actins.size(); j++){
            actins.at(j)->sound_Off_Node_Info();
        }
        cout << endl;
        cout << "***********************" << endl;
        cout << endl;
    }
    cout << "****************************************************" << endl;
    
    return;
}

//***Functions for VTK output****//
void Network::print_VTK_File(ofstream& ofs){

    ofs << "# vtk DataFile Version 3.0" << endl;
	ofs << "Point representing Actomyosin Network model" << endl;
	ofs << "ASCII" << endl << endl;
	ofs << "DATASET UNSTRUCTURED_GRID" << endl;

    //Need the total number of points/nodes for all filaments
    int num_Points = 0;

    for(unsigned int i = 0; i < filaments.size(); i++){
        num_Points += filaments.at(i)->get_Num_Actin_Nodes();
    }

    //The node/point positions
    ofs << "POINTS " << num_Points << " float64" << endl;

    vector<int> start_points;
	vector<int> end_points;
	int count = 0;

    for(unsigned int i = 0; i < filaments.size(); i++){
        start_points.push_back(count);
        filaments.at(i)->print_VTK_Points(ofs,count);
        end_points.push_back(count - 1);
    }

    ofs << endl;

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

//This function prints the node data info
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
    ofs << "NumberFilaments" << ' ' << "Ti" << endl;
    ofs << this->get_Num_Filaments() << ' ' << Ti << endl;
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

