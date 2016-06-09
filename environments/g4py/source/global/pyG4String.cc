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
// $Id: pyG4String.cc,v 1.4 2006-06-29 15:33:21 gunter Exp $
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
