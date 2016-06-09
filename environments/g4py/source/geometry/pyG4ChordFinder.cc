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
// $Id: pyG4ChordFinder.cc,v 1.4 2006-06-29 15:31:53 gunter Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ChordFinder.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4ChordFinder.hh"
#include "G4MagneticField.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4ChordFinder {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_SetVerbose, SetVerbose, 0, 1);

};

using namespace pyG4ChordFinder;

// ====================================================================
// module definition
// ====================================================================
void export_G4ChordFinder()
{
  class_<G4ChordFinder, G4ChordFinder*, boost::noncopyable>
    ("G4ChordFinder", "chord finder class", no_init)
    // constructor
    .def(init<G4MagInt_Driver*>())
    .def(init<G4MagneticField*>())
    .def(init<G4MagneticField*, G4double>())
    .def(init<G4MagneticField*, G4double, G4MagIntegratorStepper*>())
    // ---
    .def("GetDeltaChord",   &G4ChordFinder::GetDeltaChord)
    .def("SetDeltaChord",   &G4ChordFinder::SetDeltaChord)
    // ---
    .def("PrintStatistics", &G4ChordFinder::PrintStatistics)
    .def("SetVerbose",      &G4ChordFinder::SetVerbose, f_SetVerbose())
    ;
}
