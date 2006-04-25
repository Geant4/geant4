// $Id: pyG4VSensitiveDetector.cc,v 1.2 2006-04-25 08:09:44 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VSensitiveDetector.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VSensitiveDetector.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4VSensitiveDetector {

class  CB_G4VSensitiveDetector :
    public G4VSensitiveDetector,
    public wrapper<G4VSensitiveDetector> {

public:
  CB_G4VSensitiveDetector() : G4VSensitiveDetector("") { }
  CB_G4VSensitiveDetector(const G4String& name) 
    : G4VSensitiveDetector(name) { }
  ~CB_G4VSensitiveDetector() { }

  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist) {
    return get_override("ProcessHits")(&aStep, &ROhist);
  }
};

};

using namespace pyG4VSensitiveDetector;

// ====================================================================
// module definition
// ====================================================================
void export_G4VSensitiveDetector()
{
  class_<CB_G4VSensitiveDetector, boost::noncopyable>
    ("G4VSensitiveDetector", "base class of senstive detector")
    // ---
    .def(init<const G4String&>())
    // ---
    .def("Initialize",      &G4VSensitiveDetector::Initialize)
    .def("EndOfEvent",      &G4VSensitiveDetector::EndOfEvent)
    .def("clear",           &G4VSensitiveDetector::clear)
    .def("DrawAll",         &G4VSensitiveDetector::DrawAll)
    .def("PrintAll",        &G4VSensitiveDetector::PrintAll)
    .def("Hit",             &G4VSensitiveDetector::Hit)
    .def("ProcessHits", pure_virtual(&CB_G4VSensitiveDetector::ProcessHits))
    // ---
    .def("SetROgeometry",   &G4VSensitiveDetector::SetROgeometry)
    .def("GetNumberOfCollections", 
	 &G4VSensitiveDetector::GetNumberOfCollections)
    .def("GetCollectionName", &G4VSensitiveDetector::GetCollectionName)
    .def("SetVerboseLevel", &G4VSensitiveDetector::SetVerboseLevel)
    .def("Activate",        &G4VSensitiveDetector::Activate)
    .def("isActive",        &G4VSensitiveDetector::isActive)
    .def("GetName",         &G4VSensitiveDetector::GetName)
    .def("GetPathName",     &G4VSensitiveDetector::GetPathName)
    .def("GetFullPathName", &G4VSensitiveDetector::GetFullPathName)
    .def("GetROgeometry",   &G4VSensitiveDetector::GetROgeometry,
	 return_internal_reference<>())
    ;
}
