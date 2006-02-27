// $Id: pyG4VPhysicalVolume.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VPhysicalVolume.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VPVParameterisation.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4VPhysicalVolume {

const G4RotationMatrix*(G4VPhysicalVolume::*f1_GetRotation)() const
  = &G4VPhysicalVolume::GetRotation;

G4RotationMatrix*(G4VPhysicalVolume::*f2_GetRotation)()
  = &G4VPhysicalVolume::GetRotation;

};

using namespace pyG4VPhysicalVolume;

// ====================================================================
// module definition
// ====================================================================
void export_G4VPhysicalVolume()
{
  class_<G4VPhysicalVolume, G4VPhysicalVolume*, boost::noncopyable>
    ("G4VPhysicalVolume", "physical volume class", no_init)
    // ---
    .def("SetTranslation",       &G4VPhysicalVolume::SetTranslation)
    .def("GetTranslation",       &G4VPhysicalVolume::GetTranslation,
	 return_internal_reference<>())
    .def("GetObjectTranslation", &G4VPhysicalVolume::GetObjectTranslation)
    .def("GetFrameTranslation",  &G4VPhysicalVolume::GetObjectTranslation)
    // ---
    .def("SetRotation",          &G4VPhysicalVolume::SetRotation)
    .def("GetRotation",          f1_GetRotation,
      return_internal_reference<>())
    .def("GetRotation",          f2_GetRotation,
      return_internal_reference<>())
    .def("GetObjectRotation",    &G4VPhysicalVolume::GetObjectRotation,
	 return_value_policy<reference_existing_object>())
    .def("GetObjectRotationValue", &G4VPhysicalVolume::GetObjectRotationValue)
    .def("GetFrameRotation",       &G4VPhysicalVolume::GetFrameRotation,
	 return_internal_reference<>())
    // ---
    .def("SetLogicalVolume",     &G4VPhysicalVolume::SetLogicalVolume)
    .def("SetMotherLogical",     &G4VPhysicalVolume::SetMotherLogical)
    .def("GetLogicalVolume",     &G4VPhysicalVolume::GetLogicalVolume,
	 return_internal_reference<>())
    .def("GetMotherLogical",     &G4VPhysicalVolume::GetMotherLogical,
	 return_internal_reference<>())
    // ---
    .def("SetName",             &G4VPhysicalVolume::SetName)
    .def("GetName",             &G4VPhysicalVolume::GetName)
    .def("SetCopyNo",           &G4VPhysicalVolume::SetCopyNo)
    .def("GetCopyNo",           &G4VPhysicalVolume::GetCopyNo)
    // ---
    .def("IsMany",              &G4VPhysicalVolume::IsMany)
    .def("IsReplicated",        &G4VPhysicalVolume::IsReplicated)
    .def("IsParameterised",     &G4VPhysicalVolume::IsParameterised)
    .def("GetMultiplicity",     &G4VPhysicalVolume::GetMultiplicity)
    .def("GetParameterisation", &G4VPhysicalVolume::GetParameterisation,
	 return_value_policy<reference_existing_object>())
    ;
}

