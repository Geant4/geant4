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
// $Id: pyG4NistManager.cc,v 1.6 2007-11-07 08:49:31 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4NistManager.cc
//
//                                         2005 Q
// ===================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#if G4VERSION_NUMBER >= 710
#include "G4NistManager.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4NistManager {

// FindOrBuildElement
G4Element*(G4NistManager::*f1_FindOrBuildElement)(G4int, G4bool)
  = &G4NistManager::FindOrBuildElement;

G4Element*(G4NistManager::*f2_FindOrBuildElement)(const G4String&, G4bool)
  = &G4NistManager::FindOrBuildElement;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_FindOrBuildElement, 
				       FindOrBuildElement, 1, 2);

// PrintElement
void(G4NistManager::*f1_PrintElement)(const G4String&)
  = &G4NistManager::PrintElement;
void(G4NistManager::*f2_PrintElement)(G4int)
  = &G4NistManager::PrintElement;

// FindOrBuildMaterial
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_FindOrBuildMaterial, 
				       FindOrBuildMaterial, 1, 2);

// ConstructNewMaterial
#if G4VERSION_NUMBER >= 910
G4Material*(G4NistManager::*f1_ConstructNewMaterial)
  (const G4String&, const std::vector<G4String>&, 
   const std::vector<G4int>&, G4double, G4bool,
   G4State, G4double, G4double)
  = &G4NistManager::ConstructNewMaterial;

G4Material*(G4NistManager::*f2_ConstructNewMaterial)
  (const G4String&, const std::vector<G4String>&, 
   const std::vector<G4double>&, G4double, G4bool,
   G4State, G4double, G4double)
  = &G4NistManager::ConstructNewMaterial;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_ConstructNewMaterial,
				       ConstructNewMaterial, 4, 8);

#else
G4Material*(G4NistManager::*f1_ConstructNewMaterial)
  (const G4String&, const std::vector<G4String>&, 
   const std::vector<G4int>&, G4double, G4bool)
  = &G4NistManager::ConstructNewMaterial;

G4Material*(G4NistManager::*f2_ConstructNewMaterial)
  (const G4String&, const std::vector<G4String>&, 
   const std::vector<G4double>&, G4double, G4bool)
  = &G4NistManager::ConstructNewMaterial;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_ConstructNewMaterial,
				       ConstructNewMaterial, 4, 5);

#endif

#if G4VERSION_NUMBER >= 910
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_ConstructNewGasMaterial,
				       ConstructNewGasMaterial, 4, 5);
#endif

};

using namespace pyG4NistManager;

#endif

// ====================================================================
// module definition
// ====================================================================
void export_G4NistManager()
{
#if G4VERSION_NUMBER >= 710

  class_<G4NistManager, boost::noncopyable>
    ("G4NistManager", "manager class for NIST materials", no_init)
    // ---
    .def("Instance", &G4NistManager::Instance,
         return_value_policy<reference_existing_object>())
    .staticmethod("Instance")
    // ---
    .def("SetVerbose",          &G4NistManager::SetVerbose)
    .def("GetVerbose",          &G4NistManager::GetVerbose)
    // ---
#if G4VERSION_NUMBER < 910
    .def("RegisterElement",     &G4NistManager::RegisterElement)
    .def("DeRegisterElement",   &G4NistManager::DeRegisterElement)
#endif
    .def("GetElement",          &G4NistManager::GetElement,
	 return_internal_reference<>())
    .def("FindOrBuildElement",  f1_FindOrBuildElement,
	 f_FindOrBuildElement()
         [return_value_policy<reference_existing_object>()])
    .def("FindOrBuildElement",  f2_FindOrBuildElement,
	 f_FindOrBuildElement()
         [return_value_policy<reference_existing_object>()])
    .def("GetNumberOfElements", &G4NistManager::GetNumberOfElements)
    .def("GetZ",                &G4NistManager::GetZ)
    .def("GetIsotopeMass",      &G4NistManager::GetIsotopeMass)
    .def("PrintElement",        f1_PrintElement)
    .def("PrintElement",        f2_PrintElement)
    .def("PrintG4Element",      &G4NistManager::PrintG4Element)
    // ---
#if G4VERSION_NUMBER < 910
    .def("RegisterMaterial",    &G4NistManager::RegisterMaterial)
    .def("DeRegisterMaterial",  &G4NistManager::DeRegisterMaterial)
#endif
    .def("GetMaterial",         &G4NistManager::GetMaterial,
         return_value_policy<reference_existing_object>())
    .def("FindOrBuildMaterial", &G4NistManager::FindOrBuildMaterial,
	 f_FindOrBuildMaterial()
         [return_value_policy<reference_existing_object>()])
    .def("ConstructNewMaterial", f1_ConstructNewMaterial,
	 f_ConstructNewMaterial()
         [return_value_policy<reference_existing_object>()])
    .def("ConstructNewMaterial", f2_ConstructNewMaterial,
	 f_ConstructNewMaterial()
         [return_value_policy<reference_existing_object>()])
#if G4VERSION_NUMBER >= 910
    .def("ConstructNewGasMaterial", &G4NistManager::ConstructNewGasMaterial,
	 f_ConstructNewGasMaterial()
         [return_value_policy<reference_existing_object>()])
#endif
    .def("GetNumberOfMaterials", &G4NistManager::GetNumberOfMaterials)
    .def("ListMaterials",        &G4NistManager::ListMaterials)
    .def("PrintG4Material",      &G4NistManager::PrintG4Material)
    ;
#endif
}

