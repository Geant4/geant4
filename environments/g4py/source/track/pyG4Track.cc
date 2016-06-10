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
// $Id: pyG4Track.cc 66892 2013-01-17 10:57:59Z gunter $
// ====================================================================
//   pyG4Track.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Track.hh"
#include "G4VProcess.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Track()
{
  class_<G4Track, G4Track*>("G4Track", "track class")
    // ---
    .def("GetTrackID",             &G4Track::GetTrackID)
    .def("GetParentID",            &G4Track::GetParentID)
    .def("GetDynamicParticle",   &G4Track::GetDynamicParticle,
         return_internal_reference<>())
    .def("GetDefinition",          &G4Track::GetDefinition,
         return_internal_reference<>())
    .def("GetPosition",            &G4Track::GetPosition,
         return_value_policy<return_by_value>())
    .def("GetGlobalTime",          &G4Track::GetGlobalTime)
    .def("GetLocalTime",           &G4Track::GetLocalTime) 
    .def("GetProperTime",          &G4Track::GetProperTime)
    .def("GetVolume",            &G4Track::GetVolume,
         return_value_policy<reference_existing_object>())
    .def("GetMaterial",            &G4Track::GetMaterial,
	 return_value_policy<reference_existing_object>())
    .def("GetTouchable",           &G4Track::GetTouchable,
	 return_value_policy<reference_existing_object>())
    .def("GetKineticEnergy",       &G4Track::GetKineticEnergy)
    .def("GetTotalEnergy",         &G4Track::GetTotalEnergy)
    .def("GetMomentumDirection",   &G4Track::GetMomentumDirection,
	 return_value_policy<return_by_value>())
    .def("GetMomentum",            &G4Track::GetMomentum,
	 return_value_policy<return_by_value>())
    .def("GetVelocity",            &G4Track::GetVelocity)
    .def("GetPolarization",        &G4Track::GetPolarization,
	 return_value_policy<return_by_value>())
    .def("GetTrackStatus",         &G4Track::GetTrackStatus)
    .def("GetTrackLength",         &G4Track::GetTrackLength)
    .def("GetStep",                &G4Track::GetStep,
    	 return_value_policy<reference_existing_object>())
    .def("GetCurrentStepNumber",   &G4Track::GetCurrentStepNumber)
    .def("GetStepLength",          &G4Track::GetStepLength)
    .def("GetVertexPosition",      &G4Track::GetVertexPosition,
	 return_value_policy<return_by_value>())
    .def("GetVertexMomentumDirection", &G4Track::GetVertexMomentumDirection,
	 return_value_policy<return_by_value>())
    .def("GetVertexKineticEnergy", &G4Track::GetVertexKineticEnergy)
    .def("GetLogicalVolumeAtVertex", &G4Track::GetLogicalVolumeAtVertex,
    	 return_value_policy<reference_existing_object>())
    .def("GetCreatorProcess",   &G4Track::GetCreatorProcess,
    	 return_value_policy<reference_existing_object>())
    .def("GetWeight",              &G4Track::GetWeight)
    .def("SetWeight",              &G4Track::SetWeight)
    ;
}
