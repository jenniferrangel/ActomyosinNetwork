//network.h
//=================
//Include Guards
#ifndef _NETWORK_H_INCLUDED_
#define _NETWORK_H_INCLUDED_
//=========================
//forward declarations

//=======================
//Include dependencies
#include <string>
#include <vector>
#include <fstream>
#include <memory>
#include <random>
#include "params.h"
#include "coord.h"
#include "node.h"
#include "filament.h"
#include "externs.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
//=========================
// Network Class Declaration

//class Network: public enable_shared_from_this<Network> {
class Network{
    private:
        vector<shared_ptr<Filament>> filaments;

        int num_filaments;
        mt19937 gen;

    public:
        //Constructor:
        //==========================
        Network(string filename, mt19937 gen);
        Network(string filename, mt19937 gen, bool predefined_nodes);        

        //Getters & setters:
        //==========================
        void get_Filaments(vector<shared_ptr<Filament>>& filaments);
        //set and get the number of actin filaments in the network:
        void update_Num_Filaments(shared_ptr<Filament>& new_Filament);  //updates # filaments but also filaments vector
        int get_Num_Filaments(){return num_filaments;}
        mt19937 get_Random_Generator(){return gen;}
        double get_Normally_Distributed_Random_Number(double mean, double stddev);

        //Functions:
        //==========================
        //for double-checking purposes
        void sound_Off_All_Node_Info();
        void sound_Off_Neighbors();

        //Calculate forces
        void calculate_New_Forces(int Ti);

        //Update node positions via Langevin equation
        void update_Positions(int Ti);

        //Printing and dataoutput functions
        //***Functions for VTK output
        void print_VTK_File(ofstream& ofs);
        void update_VTK_Indices();

        //***Printing data output
        void locations_Output(ofstream& ofs, int Ti);
        void node_Data_Output(ofstream& ofs, int Ti);
        void filament_Data_Output(ofstream& ofs, int Ti);
        void network_Data_Output(ofstream& ofs, int Ti);

        //Destructor:
        //==========================
	    ~Network();

};

//===========================
//End of file

#endif