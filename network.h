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
        //vector<shared_ptr<Filament>> filaments;

        int num_filaments;
        mt19937 gen;

    public:
        //Constructor:
        Network(string filename, mt19937 gen);

        //Getters & setters
        int get_num_filaments(){return num_filaments;}
        mt19937 get_Random_Generator(){return gen;}

        //Functions

        //Destructor
	    ~Network();

};

//===========================
//End of file

#endif