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
// $Id: pyG4String.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4String.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4String.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
std::ostream& operator<<(std::ostream& ostr, const G4String& astr)
{
  return ostr << astr.c_str();
}

// ====================================================================
// module definition
// ====================================================================
void export_G4String()
{
  class_<G4String>("G4String", "string class")
    .def(init<const G4String&>())
    .def(init<const char*>())
    .def(self_ns::str(self))
    .def(self + self)
    .def(self += self)
    .def(self += other<const char*>())
    .def(self == self)
    .def(self == other<const char*>())
    .def(self != self)
    .def(self != other<const char*>())
    ;
  
  implicitly_convertible<G4String, const char*>();
  implicitly_convertible<const char* ,G4String>();

  implicitly_convertible<G4String, std::string>();
  implicitly_convertible<std::string ,G4String>();
}
