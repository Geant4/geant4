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
// $Id: pyG4VPhysicalVolume.cc,v 1.7 2008-03-13 07:32:18 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VPhysicalVolume.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
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
#if G4VERSION_NUMBER <= 801
    .def("GetName",             &G4VPhysicalVolume::GetName)
#else
    .def("GetName",             &G4VPhysicalVolume::GetName,
         return_value_policy<return_by_value>())
#endif
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

