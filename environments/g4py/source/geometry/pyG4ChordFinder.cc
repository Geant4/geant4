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
// $Id: pyG4ChordFinder.cc,v 1.3 2006-06-04 21:34:28 kmura Exp $
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
