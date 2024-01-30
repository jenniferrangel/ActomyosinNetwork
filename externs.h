//externs.h
//==========
//Include guards
#ifndef _EXTERNS_H_INCLUDED_
#define _EXTERNS_H_INCLUDED_

//Declares external nonconstant global variables that will be defined in main.cpp

//Actin parameters:
//============================================
//parameters for linear springs:
extern double ACTIN_SPRING_EQUI_LEN;
extern double K_LINEAR_STIFF_ACTIN;
//parameters for bending springs:
extern double K_BEND_STIFF_ACTIN;
extern double THETA_EQUI_ANGLE;
extern double ACTIN_DRAG_COEFF;
extern double VISCOSITY_CYTOPLASM;
extern double ACTIN_PERSISTENCE_LEN;

//Myosin parameters:
//============================================
extern double F_MYO_PULLING;
extern double MYO_DRAG_COEFF;

//Actomyosin parameters:
//============================================
extern double ACTIN_MYO_CONNECT_RADIUS;
extern double K_LINEAR_STIFF_ACTIN_MYO_CONN;
extern double ACTIN_MYO_CONN_SPRING_EQUI_LEN;

extern double kB;
extern double TEMPERATURE;
extern double dt;

#endif