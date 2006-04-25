// $Id: pyG4FieldManager.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4FieldManager.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4ChordFinder.hh"
#include "G4MagneticField.hh"
#include "G4Track.hh"

using namespace boost::python;

// ====================================================================
// miscs
// ====================================================================
// What a hell!
inline const G4ChordFinder* G4FieldManager::GetChordFinder() const
{
  return fChordFinder;
}

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4FieldManager {

G4ChordFinder*(G4FieldManager::*f1_GetChordFinder)()
  = &G4FieldManager::GetChordFinder;
const G4ChordFinder*(G4FieldManager::*f2_GetChordFinder)() const 
  = &G4FieldManager::GetChordFinder;

};

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

