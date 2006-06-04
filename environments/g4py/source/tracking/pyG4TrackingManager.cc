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
// $Id: pyG4TrackingManager.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
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

