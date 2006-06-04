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
// $Id: test03.cc,v 1.3 2006-06-04 21:35:59 kmura Exp $
// ====================================================================
//   test03.cc
//
//   Singleton class (w/o public constructor) can be exposed?
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

#include "AClass.hh"

using namespace boost::python;

BOOST_PYTHON_MODULE(test03) {
  class_<AClass>("AClass", "Singleton")
    .def("GetPointer", &AClass::GetPointer,
	 return_value_policy<manage_new_object>())
    ;
}

