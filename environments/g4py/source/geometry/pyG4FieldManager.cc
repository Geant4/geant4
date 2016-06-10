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
// $Id: pyG4FieldManager.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4FieldManager.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4ChordFinder.hh"
#include "G4MagneticField.hh"
#include "G4Track.hh"

using namespace boost::python;

// ====================================================================
// miscs
// ====================================================================
#if G4VERSION_NUMBER <= 801
// What a hell!
inline const G4ChordFinder* G4FieldManager::GetChordFinder() const
{
  return fChordFinder;
}
#endif

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4FieldManager {

G4ChordFinder*(G4FieldManager::*f1_GetChordFinder)()
  = &G4FieldManager::GetChordFinder;
const G4ChordFinder*(G4FieldManager::*f2_GetChordFinder)() const
  = &G4FieldManager::GetChordFinder;

}

using namespace pyG4FieldManager;

// ====================================================================
// module definition
// ====================================================================
void export_G4FieldManager()
{
  class_<G4FieldManager, G4FieldManager*, boost::noncopyable>
    ("G4FieldManager", "field manager class")
    // constructors
    .def(init<>())
    .def(init<G4Field*>())
    .def(init<G4Field*, G4ChordFinder*>())
    .def(init<G4Field*, G4ChordFinder*, G4bool>())
    .def(init<G4MagneticField*>())
    // ---
    .def("SetDetectorField",      &G4FieldManager::SetDetectorField)
    .def("GetDetectorField",      &G4FieldManager::GetDetectorField,
         return_internal_reference<>())
    .def("DoesFieldExist",        &G4FieldManager::DoesFieldExist)
    .def("CreateChordFinder",     &G4FieldManager::CreateChordFinder)
    .def("SetChordFinder",        &G4FieldManager::SetChordFinder)
    .def("GetChordFinder",        f1_GetChordFinder,
	 return_internal_reference<>())
    .def("GetChordFinder",        f2_GetChordFinder,
	 return_internal_reference<>())
    .def("ConfigureForTrack",     &G4FieldManager::ConfigureForTrack)
    .def("GetDeltaIntersection",  &G4FieldManager::GetDeltaIntersection)
    .def("GetDeltaOneStep",       &G4FieldManager::GetDeltaOneStep)
    .def("SetAccuraciesWithDeltaOneStep",
	 &G4FieldManager::SetAccuraciesWithDeltaOneStep)
    .def("SetDeltaOneStep",       &G4FieldManager::SetDeltaOneStep)
    .def("SetDeltaIntersection",  &G4FieldManager::SetDeltaIntersection)
    .def("GetMinimumEpsilonStep", &G4FieldManager::GetMinimumEpsilonStep)
    .def("SetMinimumEpsilonStep", &G4FieldManager::SetMinimumEpsilonStep)
    .def("GetMaximumEpsilonStep", &G4FieldManager::GetMaximumEpsilonStep)
    .def("SetMaximumEpsilonStep", &G4FieldManager::SetMaximumEpsilonStep)
    .def("DoesFieldChangeEnergy", &G4FieldManager::DoesFieldChangeEnergy)
    .def("SetFieldChangesEnergy", &G4FieldManager::SetFieldChangesEnergy)
    ;
}

