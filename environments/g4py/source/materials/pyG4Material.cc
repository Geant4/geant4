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
// $Id: pyG4Material.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Material.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "pyG4Version.hh"
#include "G4Material.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4Material {

// AddElement
void (G4Material::*f1_AddElement)(G4Element*, G4int)
  = &G4Material::AddElement;
void (G4Material::*f2_AddElement)(G4Element*, G4double)
  = &G4Material::AddElement;

#if G4VERSION_NUMBER >= 800
BOOST_PYTHON_FUNCTION_OVERLOADS(f_GetMaterial, G4Material::GetMaterial, 1, 2);
#endif

}

// raw pointer -> Python list conversion
list f_GetFractionVector(const G4Material* material) 
{
  list fracList;
  const G4double* fracVec= material-> GetFractionVector();
  G4int nele= material-> GetNumberOfElements();
  for(G4int i=0; i<nele; i++) {
    fracList.append(fracVec[i]);
  }
  return fracList;
}

list f_GetAtomsVector(const G4Material* material)
{
  list atomsList;
  const G4int* atomsVec= material-> GetAtomsVector();
  G4int nele= material-> GetNumberOfElements();
  for(G4int i=0; i<nele; i++) {
    atomsList.append(atomsVec[i]);
  }
  return atomsList;
}

list f_GetVecNbOfAtomsPerVolume(const G4Material* material)
{
  list nbOfAtomsPerVolumeList;
  const G4double* nbOfAtomsPerVolumeVec= material-> GetVecNbOfAtomsPerVolume();
  G4int nele= material-> GetNumberOfElements();
  for(G4int i=0; i<nele; i++) {
    nbOfAtomsPerVolumeList.append(nbOfAtomsPerVolumeVec[i]);
  }
  return nbOfAtomsPerVolumeList;
}

list f_GetAtomicNumDensityVector(const G4Material* material)
{
  list atomicNumDensityList;
  const G4double* atomicNumDensityVec= material-> GetAtomicNumDensityVector();
  G4int nele= material-> GetNumberOfElements();
  for(G4int i=0; i<nele; i++) {
    atomicNumDensityList.append(atomicNumDensityVec[i]);
  }
  return atomicNumDensityList;
}

// copy constructor is private, so ...
void Print(G4Material& mat)
{
  G4cout << mat;
}

using namespace pyG4Material;

// ====================================================================
// module definition
// ====================================================================
void export_G4Material()
{
  class_<G4Material, G4Material*, boost::noncopyable>
    ("G4Material", "material class", no_init)
    .def(init<const G4String&, G4double, G4double, G4double>())
    .def(init<const G4String&, G4double, G4int>())
    // ---
    .def("AddElement",          f1_AddElement)
    .def("AddElement",          f2_AddElement)
    .def("AddMaterial",         &G4Material::AddMaterial)
    .def("GetName",             &G4Material::GetName)
    .def("SetChemicalFormula",  &G4Material::SetChemicalFormula)
    .def("GetChemicalFormula",  &G4Material::GetChemicalFormula)
    .def("GetDensity",          &G4Material::GetDensity)
    .def("GetState",            &G4Material::GetState)
    .def("GetTemperature",      &G4Material::GetTemperature)
    .def("GetPressure",         &G4Material::GetPressure)
    // ---
    .def("GetElementVector",    &G4Material::GetElementVector,
	 return_internal_reference<>())
    .def("GetElement",          &G4Material::GetElement,
         return_value_policy<reference_existing_object>())
    .def("GetTotNbOfAtomsPerVolume",  &G4Material::GetTotNbOfAtomsPerVolume)
    .def("GetTotNbOfElectPerVolume",  &G4Material::GetTotNbOfElectPerVolume)
    .def("GetFractionVector",         f_GetFractionVector)
    .def("GetAtomsVector",            f_GetAtomsVector)
    .def("GetVecNbOfAtomsPerVolume",  f_GetVecNbOfAtomsPerVolume)
    .def("GetAtomicNumDensityVector", f_GetAtomicNumDensityVector)
    // ----
    .def("GetElectronDensity",        &G4Material::GetElectronDensity)
    .def("GetRadlen",                 &G4Material::GetRadlen)
    .def("GetNuclearInterLength",     &G4Material::GetNuclearInterLength)
    .def("GetIonisation",             &G4Material::GetIonisation,
	 return_internal_reference<>())
    .def("GetSandiaTable",            &G4Material::GetSandiaTable,
	 return_internal_reference<>())
    // ---
    .def("GetZ",                  &G4Material::GetZ)
    .def("GetA",                  &G4Material::GetA)
    .def("SetMaterialPropertiesTable", &G4Material::SetMaterialPropertiesTable)
    .def("GetMaterialPropertiesTable", &G4Material::GetMaterialPropertiesTable,
	 return_internal_reference<>())
    .def("GetMaterialTable",      &G4Material::GetMaterialTable,
	 return_value_policy<reference_existing_object>())
    .staticmethod("GetMaterialTable")
    .def("GetNumberOfMaterials",  &G4Material::GetNumberOfMaterials)
    .staticmethod("GetNumberOfMaterials")
    .def("GetIndex",              &G4Material::GetIndex)
#if G4VERSION_NUMBER >= 800
    .def("GetMaterial",           &G4Material::GetMaterial,
	 f_GetMaterial()
	 [return_value_policy<reference_existing_object>()])
#else
    .def("GetMaterial",           &G4Material::GetMaterial,
         return_value_policy<reference_existing_object>())
#endif
    .staticmethod("GetMaterial")
    // ---
    //.def(self_ns::str(self))
    .def("Print", Print)
    ;

  // ---
  enum_<G4State>("G4State")
    .value("kStateUndefined", kStateUndefined)
    .value("kStateSolid",     kStateSolid)
    .value("kStateLiquid",    kStateLiquid)
    .value("kStateGas",       kStateGas)
    ;
}

