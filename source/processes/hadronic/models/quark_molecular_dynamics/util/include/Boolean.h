/*
Example programs from the book Scientific and Engineering Programming
in C++: An Introduction with Advanced Techniques and Examples,
Addison-Wesley, 1994.
 
                 (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
                 CORPORATION 1994.  ALL RIGHTS RESERVED.

See README file for further details.
*/
#ifndef BOOLEANH
#define BOOLEANH

#include "globals.hh"
#include "g4std/iostream"

class Boolean {
  // I/O

friend G4std::ostream& operator<<(G4std::ostream& s, Boolean  b);
friend G4std::istream& operator>>(G4std::istream& s, Boolean& b);

public:
    // Constants
    enum constants { False = 0, True = 1 };

    // Construction.

    Boolean()                      {}   // Construct uninitialized.
    Boolean(int i) :    v(i != 0)  {}   // Construct and initialize to (i != 0).
    Boolean(float f) :  v(f != 0)  {}   // Construct and initialize to (f != 0).
    Boolean(double d) : v(d != 0)  {}   // Construct and initialize to (d != 0).
    Boolean(void* p) :  v(p != 0)  {}   // Construct and initialize to (p != 0).

    // Conversion.
    operator int() const{ return v; }   // To allow "if (boolean-value)..."

    // Negation.
    Boolean operator!() const { return !v; }


private:
    char v;
};


#endif
