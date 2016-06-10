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
// $Id: test08.cc 86624 2014-11-14 08:59:40Z kmura $
// ====================================================================
//   test08.cc
//
//   test for static member function
//
//                                         2005 Q
// ====================================================================

class AClass {
public:
  AClass() { }
  ~AClass() { }
  static int AMethod() { return 1; } 

};

// Boost.Python...
#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(test08)
{
  class_<AClass>( "AClass", "a class")
    .def(init<>())
    .def("AMethod", &AClass::AMethod)
    .staticmethod("AMethod")
    ;
}

