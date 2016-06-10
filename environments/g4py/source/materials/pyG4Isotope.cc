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
// $Id: pyG4Isotope.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4Isotope.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#include "pyG4indexing.hh"
#include "G4Isotope.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4Isotope {

BOOST_PYTHON_FUNCTION_OVERLOADS(f_GetIsotope, G4Isotope::GetIsotope, 1, 2)

// copy constructor is private, so ...
void Print(G4Isotope& iso)
{
  G4cout << iso;
}

}

using namespace pyG4Isotope;

// ====================================================================
// module definition
// ====================================================================
void export_G4Isotope()
{
  class_<G4Isotope, G4Isotope*, boost::noncopyable>
    ("G4Isotope", "isotope class", no_init)
    // constructors
    .def(init<const G4String&, G4int, G4int>())
    .def(init<const G4String&, G4int, G4int, G4double>())
    // ---
    .def("GetName",             &G4Isotope::GetName,
         return_value_policy<reference_existing_object>())
    .def("SetName",             &G4Isotope::SetName)
    .def("GetZ",                &G4Isotope::GetZ)
    .def("GetN",                &G4Isotope::GetN)
    .def("GetA",                &G4Isotope::GetA)
    .def("GetIsotope",          &G4Isotope::GetIsotope,
         f_GetIsotope()
         [return_value_policy<reference_existing_object>()])
    .staticmethod("GetIsotope")
    .def("GetIsotopeTable",     &G4Isotope::GetIsotopeTable,
         return_value_policy<reference_existing_object>())
    .staticmethod("GetIsotopeTable")
    .def("GetNumberOfIsotopes", &G4Isotope::GetNumberOfIsotopes)
    .staticmethod("GetNumberOfIsotopes")

    .def("GetIndex",            &G4Isotope::GetIndex)
    // ---
    .def("Print", Print)
    .def(self == self)
    .def(self != self)
    ;

  // G4IsotopeTable
  class_<G4IsotopeTable> ("G4IsotopeTable", "isotope table")
    .def(vector_indexing_suite<G4IsotopeTable>())
    .def(self_ns::str(self))
    ;
}
