//coord.h
#ifndef _COORD_H_INCLUDED_
#define _COORD_H_INCLUDED_
//======================
//forward dependencies

//======================
//include dependencies
#include <iostream>
using namespace std;
//======================

class Coord {
    //Atributes:
    protected:
        double x;
        double y;

    //Methods:
    public:
        //Constructors:
        Coord(); // default => sets it to (0,0)
		Coord(double x, double y);
		Coord(const Coord& c);  //copy constructor => creates a new object as a copy of an existing object of the same class; called automatically when an object is passed by value or is copied

        //Getters & setters:
        double get_X() const;
		double get_Y() const;

        //Operator overloading:
        void operator=(const Coord& c);
		Coord operator-(const Coord& c) const;
		Coord operator+(const Coord& c) const;
		void operator+=(const Coord& c);
		void operator-=(const Coord& c);
		Coord operator/(const double d) const;
		Coord operator*(const double d) const;
		bool operator==(const Coord& c); 
		bool operator!=(const Coord&c);

        //Higher math fuctions:
        double dot(const Coord& c) const;
		double cross(const Coord& c) const;
        double length() const;
        Coord distribute(const Coord& c) const;
        Coord projectOnto(const Coord& c) const;
		Coord perpVector() const;

        //Display functions:
        friend ostream& operator<<(ostream& os, const Coord& c);

};

#endif