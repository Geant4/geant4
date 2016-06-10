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
// $Id: test02.cc 86749 2014-11-17 15:03:05Z gcosmo $
// ====================================================================
//   test02.cc
//
//   test of simple inheritance
//   - pure virtual class
//   - class inheritance
//   - function override
//   - polymorphic behavior
//   - pointer argument and down cast
//   - const flags
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

#include "XBase.hh"
#include "AClass.hh"
#include "BClass.hh"

using namespace boost::python;

BOOST_PYTHON_MODULE(test02) {
  class_<XBase, boost::noncopyable>("XBase", "Base Class", no_init)
    .add_property("ival", &XBase::GetIVal, &XBase::SetIVal)
    .add_property("dval", &XBase::GetDVal, &XBase::SetDVal)
    .def("AMethod", &XBase::AMethod)
    .def("VMethod", &XBase::VMethod)
    ;
  
  class_<AClass, bases<XBase> >( "AClass", "Derived Class A")
    .def(init<>())
    .def("AMethod", &AClass::AMethod)
    ;

  class_<BClass, bases<XBase> >( "BClass", "Derived Class B")
    .def(init<>())
    .def("AMethod", &BClass::AMethod)
    ;
}

