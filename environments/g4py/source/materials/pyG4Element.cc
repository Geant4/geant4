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
// $Id: pyG4Element.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Element.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
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

};

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
    .def("GetName",             &G4Element::GetName)
    .def("GetSymbol",           &G4Element::GetSymbol)
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
    .def("GetCountUse",         &G4Element::GetCountUse)
    .def("increaseCountUse",    &G4Element::increaseCountUse)
    .def("decreaseCountUse",    &G4Element::decreaseCountUse)
    .def("GetIndexZ",           &G4Element::GetIndexZ)
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

