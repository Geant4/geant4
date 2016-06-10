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
// $Id: pyG4Element.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4Element.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#include "G4Element.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4Element {

// raw pointer -> Python list conversion
list f_GetRelativeAbundanceVector(const G4Element* element)
{
  list aList;
  const G4double* aVec= element-> GetRelativeAbundanceVector();
  G4int niso= element-> GetNumberOfIsotopes();
  for(G4int i=0; i<niso; i++) {
    aList.append(aVec[i]);
  }
  return aList;
}

// copy constructor is private, so ...
void Print(G4Element& ele)
{
  std::cout << ele;  // problem with G4cout. (delayed message)
}

}

using namespace pyG4Element;

// ====================================================================
// module definition
// ====================================================================
void export_G4Element()
{
  class_<G4Element, G4Element*, boost::noncopyable>
    ("G4Element", "element class", no_init)
    // constructors
    .def(init<const G4String&, const G4String&, G4double, G4double>())
    .def(init<const G4String&, const G4String&, G4int>())
    // ---
    .def("AddIsotope",          &G4Element::AddIsotope)
    .def("GetName",             &G4Element::GetName,
         return_value_policy<reference_existing_object>())
    .def("GetSymbol",           &G4Element::GetSymbol,
         return_value_policy<reference_existing_object>())
    .def("SetName",             &G4Element::SetName)
    .def("GetZ",                &G4Element::GetZ)
    .def("GetN",                &G4Element::GetN)
    .def("GetA",                &G4Element::GetA)
    .def("GetNbOfAtomicShells", &G4Element::GetNbOfAtomicShells)
    .def("GetAtomicShell",      &G4Element::GetAtomicShell)
    .def("GetNumberOfIsotopes", &G4Element::GetNumberOfIsotopes)
    .def("GetIsotopeVector",    &G4Element::GetIsotopeVector,
	 return_internal_reference<>())
    .def("GetRelativeAbundanceVector", f_GetRelativeAbundanceVector)
    .def("GetIsotope",          &G4Element::GetIsotope,
	 return_value_policy<reference_existing_object>())
    .def("GetElementTable",     &G4Element::GetElementTable,
	 return_value_policy<reference_existing_object>())
    .staticmethod("GetElementTable")
    .def("GetNumberOfElements", &G4Element::GetNumberOfElements)
    .staticmethod("GetNumberOfElements")
    .def("GetIndex",            &G4Element::GetIndex)
    .def("GetElement",          &G4Element::GetElement,
         return_value_policy<reference_existing_object>())
    .staticmethod("GetElement")
    .def("GetfCoulomb",         &G4Element::GetfCoulomb)
    .def("GetfRadTsai",         &G4Element::GetfRadTsai)
    .def("GetIonisation",       &G4Element::GetIonisation,
	 return_internal_reference<>())
    // ---
    .def("Print", Print)
    .def(self == self)
    .def(self != self)
    ;

}
