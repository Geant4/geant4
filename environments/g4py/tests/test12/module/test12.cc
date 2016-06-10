//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: test12.cc 86749 2014-11-17 15:03:05Z gcosmo $
// ====================================================================
//   test12.cc
//
//   test for indexing a STL vector
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#if BOOST_VERSION < 103300
#include "vector_indexing_suite.hpp"
#else 
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif

#include <vector>
#include <iostream>

class AClass {
private:
  int ival;

public:
  AClass() : ival(0) {
  }

  AClass(const AClass& right) {
    ival= right.ival;
    std::cout << "*** copy constructor is called" << std::endl;
  }

  ~AClass() { }

  void SetIVal(int i) { ival= i; }
  int GetIVal() { return ival; }

  void Print() const {
    std::cout << "*** @" << this << ": i=" 
	      << ival << std::endl;
  }

};

typedef std::vector<AClass*> AVector;


void PrintVector(const AVector& vec) {
  std::cout << "*** size of AVector= " << vec.size() << std::endl;
  for (int i=0; i< vec.size(); i++) {
    vec[i]-> Print();
  }
}


// Boost.Python...

using namespace boost::python;

BOOST_PYTHON_MODULE(test12)
{
  class_<AClass, AClass*>("AClass", "a class")
    .add_property("ival", &AClass::GetIVal, &AClass::SetIVal)
    .def("Print",         &AClass::Print)
    ;

  class_<AVector> ("AVector", "AClass vector")
    .def(vector_indexing_suite<AVector>())
    ;

  def("PrintVector", PrintVector);
}

