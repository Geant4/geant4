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
// $Id: pyG4VSensitiveDetector.cc,v 1.3 2006-06-04 21:34:28 kmura Exp $
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
