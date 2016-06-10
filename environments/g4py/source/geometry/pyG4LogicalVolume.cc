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
// $Id: pyG4LogicalVolume.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4LogicalVolume.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VSolid.hh"
#include "G4FieldManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4UserLimits.hh"
#include "G4SmartVoxelHeader.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4FastSimulationManager.hh"
#include "G4VisAttributes.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4LogicalVolume {

void(G4LogicalVolume::*f1_SetVisAttributes)(const G4VisAttributes*)
  = &G4LogicalVolume::SetVisAttributes;

void(G4LogicalVolume::*f2_SetVisAttributes)(const G4VisAttributes&)
  = &G4LogicalVolume::SetVisAttributes;

G4VSolid*(G4LogicalVolume::*f1_GetSolid)() const = &G4LogicalVolume::GetSolid;

void(G4LogicalVolume::*f1_SetSolid)(G4VSolid*)
  = &G4LogicalVolume::SetSolid;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetMass, GetMass, 0, 3)

}

using namespace pyG4LogicalVolume;

// ====================================================================
// module definition
// ====================================================================
void export_G4LogicalVolume()
{
  class_<G4LogicalVolume, G4LogicalVolume*, boost::noncopyable>
    ("G4LogicalVolume", "logical volume class", no_init)
    // constructors
    .def(init<G4VSolid*, G4Material*, const G4String& >())
    .def(init<G4VSolid*, G4Material*, const G4String&,
         G4FieldManager* >())
    .def(init<G4VSolid*, G4Material*, const G4String&,
         G4FieldManager*, G4VSensitiveDetector* >())
    .def(init<G4VSolid*, G4Material*, const G4String&,
         G4FieldManager*, G4VSensitiveDetector*,
         G4UserLimits* >())
    .def(init<G4VSolid*, G4Material*, const G4String&,
         G4FieldManager*, G4VSensitiveDetector*,
         G4UserLimits*, G4bool >())
    // ---
    .def("GetName",         &G4LogicalVolume::GetName)
    .def("SetName",         &G4LogicalVolume::SetName)
    // ---
    .def("GetNoDaughters",  &G4LogicalVolume::GetNoDaughters)
    .def("GetDaughter",     &G4LogicalVolume::GetDaughter,
	 return_internal_reference<>())
    .def("AddDaughter",     &G4LogicalVolume::AddDaughter)
    .def("IsDaughter",      &G4LogicalVolume::IsDaughter)
    .def("IsAncestor",      &G4LogicalVolume::IsAncestor)
    .def("RemoveDaughter",  &G4LogicalVolume::RemoveDaughter)
    .def("ClearDaughters",  &G4LogicalVolume::ClearDaughters)
    .def("TotalVolumeEntities", &G4LogicalVolume::TotalVolumeEntities)
    // ----
    .def("GetSolid",        f1_GetSolid,
     return_internal_reference<>())
    .def("SetSolid",        f1_SetSolid)
    .def("GetMaterial",     &G4LogicalVolume::GetMaterial,
	 return_internal_reference<>())
    .def("SetMaterial",     &G4LogicalVolume::SetMaterial)
    .def("UpdateMaterial",  &G4LogicalVolume::UpdateMaterial)
    // ---
    .def("GetMass",         &G4LogicalVolume::GetMass, f_GetMass())
    .def("GetFieldManager", &G4LogicalVolume::GetFieldManager,
	 return_internal_reference<>())
    .def("SetFieldManager", &G4LogicalVolume::SetFieldManager)
    .def("GetSensitiveDetector", &G4LogicalVolume::GetSensitiveDetector,
	 return_internal_reference<>())
    .def("GetUserLimits",   &G4LogicalVolume::GetUserLimits,
	 return_internal_reference<>())
    .def("SetUserLimits",   &G4LogicalVolume::SetUserLimits)
    // ---
    .def("GetVoxelHeader",  &G4LogicalVolume::GetVoxelHeader,
	 return_internal_reference<>())
    .def("SetVoxelHeader",  &G4LogicalVolume::SetVoxelHeader)
    .def("GetSmartless",    &G4LogicalVolume::GetSmartless)
    .def("SetSmartless",    &G4LogicalVolume::SetSmartless)
    .def("IsToOptimise",    &G4LogicalVolume::IsToOptimise)
    .def("SetOptimisation", &G4LogicalVolume::SetOptimisation)
    // ---
    .def("IsRootRegion",    &G4LogicalVolume::IsRootRegion)
    .def("SetRegionRootFlag", &G4LogicalVolume::SetRegionRootFlag)
    .def("IsRegion",        &G4LogicalVolume::IsRegion)
    .def("SetRegion",       &G4LogicalVolume::SetRegion)
    .def("GetRegion",       &G4LogicalVolume::GetRegion,
	 return_internal_reference<>())
    .def("PropagateRegion", &G4LogicalVolume::PropagateRegion)
    .def("GetMaterialCutsCouple", &G4LogicalVolume::GetMaterialCutsCouple,
	 return_internal_reference<>())
    .def("SetMaterialCutsCouple", &G4LogicalVolume::SetMaterialCutsCouple)
    // ---
    .def("GetVisAttributes", &G4LogicalVolume::GetVisAttributes,
	 return_internal_reference<>())
    .def("SetVisAttributes", f1_SetVisAttributes)
    .def("SetVisAttributes", f2_SetVisAttributes)
    // ---
    .def("GetFastSimulationManager",
	 &G4LogicalVolume::GetFastSimulationManager,
	 return_internal_reference<>())
    // ---
    .def("SetBiasWeight",  &G4LogicalVolume::SetBiasWeight)
    .def("GetBiasWeight",  &G4LogicalVolume::GetBiasWeight)
    ;
}
