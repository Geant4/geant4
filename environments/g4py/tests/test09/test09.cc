//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: test09.cc,v 1.3 2006-06-04 21:36:00 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   test09.cc
//
//   test for operators
//
//                                         2005 Q
// ====================================================================
#include <iostream>

class AClass {
private:
  int ival;

public:
  AClass() :ival(999) { }
  AClass(int i) :ival(i) { }
  ~AClass() { }

  void SetIVal(int i) { ival= i; }
  int GetIVal() const { return ival; }

  AClass operator+(const AClass& aclass) {
    AClass atemp;
    atemp.ival= ival + aclass.ival;    
    return atemp;
  }

  AClass& operator+=(const AClass& aclass) {
    ival+= aclass.ival;
    return *this;
  }

  bool operator==(const AClass& aclass) const {
    if(ival == aclass.ival) return true;
    return false;
  }
};

std::ostream& operator<<(std::ostream& ostr, const AClass& aclass)
{
  return ostr << aclass.GetIVal();
}


// Boost.Python...
#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(test09)
{
  class_<AClass>( "AClass", "a class")
    .def(init<>())
    .def(init<int>())
    .add_property("ival", &AClass::GetIVal, &AClass::SetIVal)
    .def(self + self)
    .def(self += self)
    .def(self == self)
    .def(self_ns::str(self))
    ;
}

