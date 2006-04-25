// $Id: pyG4TrackingManager.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4TrackingManager.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4TrackingManager.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4TrackingManager()
{
  class_<G4TrackingManager>
    ("G4TrackingManager", "tracking manager class")
    // ---
    // Note that exposed items are limited, because this class object
    // is mainly for internal uses.
    .def("GetStoreTrajectory", &G4TrackingManager::GetStoreTrajectory)
    .def("SetStoreTrajectory", &G4TrackingManager::SetStoreTrajectory)
    .def("SetVerboseLevel",    &G4TrackingManager::SetVerboseLevel)
    .def("GetVerboseLevel",    &G4TrackingManager::GetVerboseLevel)
    .def("SetUserTrackInformation", 
	 &G4TrackingManager::SetUserTrackInformation)
    ;
}

