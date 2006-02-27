// $Id: pyG4Track.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Track.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Track.hh"

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
    //.def("GetVolume",            &G4Track::GetVolume,
    //return_value_policy<reference_existing_object>())
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
    .def("GetCurrentStepNumber",   &G4Track::GetCurrentStepNumber)
    .def("GetStepLength",          &G4Track::GetStepLength)
    .def("GetVertexPosition",      &G4Track::GetVertexPosition,
	 return_value_policy<return_by_value>())
    .def("GetVertexMomentumDirection", &G4Track::GetVertexMomentumDirection,
	 return_value_policy<return_by_value>())
    .def("GetVertexKineticEnergy", &G4Track::GetVertexKineticEnergy)
    //.def("GetLogicalVolumeAtVertex", &G4Track::GetLogicalVolumeAtVertex,
    //	 return_value_policy<reference_existing_object>())
    //.def("GetCreatorProcess",   &G4Track::GetCreatorProcess,
    //	 return_value_policy<reference_existing_object>())
    .def("GetWeight",              &G4Track::GetWeight)
    .def("SetWeight",              &G4Track::SetWeight)
    ;
}
