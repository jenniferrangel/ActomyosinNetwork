//coord.cpp
//=========================
//forward dependencies

//=========================
//include dependencies
#include <iostream>
#include <omp.h>
#include <unistd.h>
#include <math.h>
#include "coord.h"
//=========================
using namespace std;

// Coord Member Functions

// Constructors:
Coord::Coord() {
	x = 0;
	y = 0;
}

Coord::Coord(double x, double y) {
	this->x = x;
	this->y = y;
}

Coord::Coord(const Coord& c) {
	x = c.get_X();
	y = c.get_Y();
}

//Getters & setters:
double Coord::get_X() const{return x;}
double Coord::get_Y() const{return y;}

//Operator overloading:
void Coord::operator=(const Coord& c){
    x = c.get_X();
    y = c.get_Y();
}

Coord Coord::operator-(const Coord& c) const{
    Coord q( x - c.get_X(), y - c.get_Y() );
	return q;
}

Coord Coord::operator+(const Coord& c) const{
    Coord q( x + c.get_X(), y + c.get_Y() );
	return q;
}

void Coord::operator+=(const Coord& c){
    x += c.get_X();
	y += c.get_Y();
}

void Coord::operator-=(const Coord& c){
    x -= c.get_X();
	y -= c.get_Y();
}

Coord Coord::operator/(const double d) const{
    Coord q( x / d, y / d);
	return q;
}

Coord Coord::operator*(const double d) const{
    Coord q( x * d, y * d);
	return q;
}

bool Coord::operator==(const Coord& c){
    if((x == c.get_X())&&(y==c.get_Y())){
		return true;
	}
	else {
		return false;
	}
}

bool Coord::operator!=(const Coord&c){
    if((x != c.get_X()) || (y != c.get_Y())){
		return true;
	}
	else {
		return false;
	}
}

//Higher math fuctions:
double Coord::dot(const Coord& c) const{
    return ( (x * c.get_X()) + (y * c.get_Y()) );  //a.b = sum(a_i.b_i) = a1b1 + a2b2 + ... + anbn
}

double Coord::cross(const Coord& c) const{
    return (x * c.get_Y()) -  (y * c.get_X());
}

double Coord::length() const{
    return sqrt( (x * x) + (y * y) );
}

Coord Coord::distribute(const Coord& c) const{
    Coord q(x*c.get_X(),y*c.get_Y());
	return q;
}

//Returns the projection of this onto c.
Coord Coord::projectOnto(const Coord& c) const { 
	Coord proj;
	if ( c.length() == 0 ) { 
		Coord temp(0,0);
		proj = temp;
	} else { 
		proj =  c * (c.dot(*this) / pow(c.length(),2) );
	}
	return proj;
}

//Returns a unit-length perpendicular vector to this with nonnegative y component.
Coord Coord::perpVector() const { 
	if (x == 0) {
		if (y == 0) {
			Coord q(0,0);
			return q;
		} else {
			Coord q(1,0);
			return q;
		}
	} else {
		Coord q(-y/x,1);
		return q/q.length();
	}

}
 
//Display functions:
ostream& operator<<(ostream& os, const Coord& c) {
	os << '(' << c.get_X() << ',' << c.get_Y() << ')';
	return os;
}