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
// $Id: test05.cc,v 1.4 2006-06-29 15:38:41 gunter Exp $
// ====================================================================
//   test05.cc
//
//   test of function overloading
//   - constructor
//   - default arguments
//   - auto-overloading
//   - docstrings
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

#include "XFunc.hh"
#include "AClass.hh"

using namespace boost::python;

BOOST_PYTHON_FUNCTION_OVERLOADS(f_func2, func2, 1, 2);
//int(*func2_1)(int) = func2;
//int(*func2_2)(int,int) = func2;

int(*func3_1)(int)= func3;
int(*func3_2)(double)= func3;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_AMethod, AMethod, 0, 2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_CMethod, CMethod, 1, 3);

int(AClass::*f1_BMethod)() = &AClass::BMethod;
double(AClass::*f2_BMethod)(double) = &AClass::BMethod;


BOOST_PYTHON_MODULE(test05) {
  // this is a temporal solution
  def("func2", (int(*)(int,int))0, f_func2());
  //def("func2", func2_1);
  //def("func2", func2_2);

  def("func3", func3_2);
  def("func3", func3_1);

  class_<AClass>( "AClass", "a class")
    .def(init<>())
    .def(init<int, optional<double> >())
    .add_property("i", &AClass::GetIVal, &AClass::SetIVal)
    .def("AMethod", (int(AClass::*)(int,double))0,
	 f_AMethod(args("i","d"), "amethod"))
    .def("BMethod", f1_BMethod)
    .def("BMethod", f2_BMethod, args("i","d"), "bmethod2")
    .def("CMethod", &AClass::CMethod,
	 f_CMethod((arg("i"), arg("d1")=1., arg("d2")=2.),
		   "cmethod: return i*d1*d2")
      )    
    ;
}

