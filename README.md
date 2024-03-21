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