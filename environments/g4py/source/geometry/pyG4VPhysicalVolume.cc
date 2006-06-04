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
// $Id: pyG4VPhysicalVolume.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
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

