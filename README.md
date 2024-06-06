# ActomyosinNetwork

This repository contains C++ code for modeling an actomysoin network composed of actin filaments and nonmuscle myosin II minifilaments.

Software environment required:

    *CMAKE ------ Build the system.
    *Paraview --- (Optional) Visualization software of vtk files for animation purposes.

***********************************************************************

Parameter files:

    * all parameters can be found in main.cpp, externs.h, params.h
    * params.h: contains the parameters that are calibrated, known or will not be changed (constantly or at all)
    * main.cpp & externs.h: parameters in main.cpp line up with the definitions in externs.h
        - the values of parameters are defined in main
        - externs.h makes the variable an external variable and allows us to use it outside of main

***********************************************************************

Initial conditions (IC):

    (1) txt file that only includes the position of the first node but the code generates the rest of the nodes
        - Ex: actomyo_initial_condition.txt

    (2) txt file with positions of all nodes in the filament predefined
    params.h: contains the parameters that are calibrated, known or will not be changed (constantly or at all)
        - Ex: actomyo_nodes_initial.txt
    
    Instructions: In main.cpp
        - set PREDEFINED_INITIAL_NODES = false if using IC (1) 
        - set PREDEFINED_INITIAL_NODES = true if using IC (2) 
        - uncomment/comment respective lines under the section "//instantiate the network:"

    Can test different network configurations: Open folder "Make_Actomyosin_Network_Scripts"
        - This folder contains 3 different subfolders
            * Network_Straight_Filaments
            * Network_Random_Alignment
            * Network_Vertical_Horizontal_Random_Alignment
        - Each subfolder includes
            *sample txt file with a network
            *matlab code to generate each type of network 

***********************************************************************

To compile and run the code use the following commands:

    (1) make clean
    (2) make all
    (3) ./program 

    Note: after compiling 2 folders will be created (Animation and DataOutput)
        *inside DataOutput create subfolders:
            -Locations
            -Node_Data
            -Filament_Data
            -Tissue_Data

        *Then run using ./program Animation DataOutput/Locations DataOutput Node_Data DataOutput/Filament_Data DataOutput/Tissue_Data


***********************************************************************