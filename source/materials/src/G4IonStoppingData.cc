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
//
// ===========================================================================
// GEANT4 class source file
//
// Class:                G4IonStoppingData
//
// Base class:           G4VIonDEDXTable 
// 
// Author:               Anton Lechner (Anton.Lechner@cern.ch)
//
// First implementation: 03. 11. 2009
//
// Modifications:
// 25.10.2010 V.Ivanchenko fixed warnings reported by the Coverity tool
// 25.10.2011: new scheme for G4Exception  (mma)
//
//
// Class description: Class which can read ion stopping power data from
//                    $G4LEDATA/ion_stopping_data
//
// Comments:
//
// =========================================================================== 
//

#include <fstream>
#include <sstream>
#include <iomanip>

#include "G4IonStoppingData.hh" 
#include "G4PhysicsVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// #########################################################################

G4IonStoppingData::G4IonStoppingData(const G4String& dir, G4bool val) :
  subDir(dir), fICRU90(val) {
}

// #########################################################################

G4IonStoppingData::~G4IonStoppingData() {
  ClearTable();
}

// #########################################################################

G4bool G4IonStoppingData::IsApplicable(
         G4int atomicNumberIon,  // Atomic number of ion
         G4int atomicNumberElem  // Atomic number of elemental material
				    )
{
  G4IonDEDXKeyElem key = std::make_pair(atomicNumberIon, atomicNumberElem);

  auto iter = dedxMapElements.find(key);

  return iter != dedxMapElements.end();
}

// #########################################################################

G4bool G4IonStoppingData::IsApplicable(
         G4int atomicNumberIon,         // Atomic number of ion
         const G4String& matIdentifier  // Name or chemical formula of material
				    )
{
  G4IonDEDXKeyMat key = std::make_pair(atomicNumberIon, matIdentifier);

  auto iter = dedxMapMaterials.find(key);

  return iter != dedxMapMaterials.end();
}

// #########################################################################

G4PhysicsVector* G4IonStoppingData::GetPhysicsVector(
         G4int atomicNumberIon,        // Atomic number of ion
         G4int atomicNumberElem        // Atomic number of elemental material
				    )
{
  G4IonDEDXKeyElem key = std::make_pair(atomicNumberIon, atomicNumberElem);

  auto iter = dedxMapElements.find(key);

  return (iter != dedxMapElements.end()) ? iter->second : nullptr; 
}

// #########################################################################

G4PhysicsVector*  G4IonStoppingData::GetPhysicsVector(
         G4int atomicNumberIon,        // Atomic number of ion
         const G4String& matIdentifier // Name or chemical formula of material
				    )
{
  G4IonDEDXKeyMat key = std::make_pair(atomicNumberIon, matIdentifier);

  auto iter = dedxMapMaterials.find(key);

  return (iter != dedxMapMaterials.end()) ? iter->second : nullptr; 
}

// #########################################################################

G4double G4IonStoppingData::GetDEDX(
         G4double kinEnergyPerNucleon, // Kinetic energy per nucleon
         G4int atomicNumberIon,        // Atomic number of ion
         G4int atomicNumberElem        // Atomic number of elemental material
				  )
{
  G4IonDEDXKeyElem key = std::make_pair(atomicNumberIon, atomicNumberElem);

  auto iter = dedxMapElements.find(key);

  return ( iter != dedxMapElements.end()) ?
    (iter->second)->Value( kinEnergyPerNucleon) : 0.0;
}

// #########################################################################

G4double G4IonStoppingData::GetDEDX(
         G4double kinEnergyPerNucleon, // Kinetic energy per nucleon
         G4int atomicNumberIon,        // Atomic number of ion
         const G4String& matIdentifier // Name or chemical formula of material
				  )
{
  G4IonDEDXKeyMat key = std::make_pair(atomicNumberIon, matIdentifier);

  auto iter = dedxMapMaterials.find(key);

  return (iter != dedxMapMaterials.end()) ?
    (iter->second)->Value(kinEnergyPerNucleon) : 0.0;
}

// #########################################################################

G4bool G4IonStoppingData::AddPhysicsVector(
        G4PhysicsVector* physicsVector, // Physics vector
	G4int atomicNumberIon,          // Atomic number of ion
        const G4String& matIdentifier   // Name of elemental material
				      ) 
{
  if(physicsVector == nullptr) {
    G4Exception ("G4IonStoppingData::AddPhysicsVector() for material", 
		 "mat037", FatalException, 
		 "Pointer to vector is null-pointer.");
    return false;
  }

  if(matIdentifier.empty()) {
    G4Exception ("G4IonStoppingData::AddPhysicsVector() for material", 
                 "mat038", FatalException, "Invalid name of the material.");
    return false;
  }

  if(atomicNumberIon <= 0) {
    G4Exception ("G4IonStoppingData::AddPhysicsVector() for material", 
                 "mat039", FatalException, "Illegal atomic number.");
    return false;
  }

  G4IonDEDXKeyMat mkey = std::make_pair(atomicNumberIon, matIdentifier);

  if(dedxMapMaterials.count(mkey) == 1) {
    G4ExceptionDescription ed;
    ed << "Vector with Z1 = " << atomicNumberIon << ", mat = " 
       << matIdentifier
       << "already exists. Remove first before replacing.";
    G4Exception ("G4IonStoppingData::AddPhysicsVector() for material", 
                 "mat040", FatalException, ed);
    return false;
  }

  dedxMapMaterials[mkey] = physicsVector;

  return true;
}

// #########################################################################

G4bool G4IonStoppingData::AddPhysicsVector(
        G4PhysicsVector* physicsVector, // Physics vector
	G4int atomicNumberIon,          // Atomic number of ion
        G4int atomicNumberElem          // Atomic number of elemental material
				      ) 
{
  if(physicsVector == nullptr) {
    G4Exception ("G4IonStoppingData::AddPhysicsVector() for element", "mat037", 
		 FatalException, "Pointer to vector is null-pointer.");
     return false;
  }

  if(atomicNumberIon <= 0) {
    G4Exception ("G4IonStoppingData::AddPhysicsVector() for element", "mat038", 
		 FatalException, "Invalid ion number.");
    return false;
  }

  if(atomicNumberElem <= 0) {
    G4Exception ("G4IonStoppingData::AddPhysicsVector() for element", "mat039", 
		 FatalException, "Illegal atomic number.");
    return false;
  }

  G4IonDEDXKeyElem key = std::make_pair(atomicNumberIon, atomicNumberElem);

  if(dedxMapElements.count(key) == 1) {
    G4ExceptionDescription ed;
    ed << "Vector with Z1 = " << atomicNumberIon << ", Z= " 
       << atomicNumberElem
       << "already exists. Remove first before replacing.";
    G4Exception ("G4IonStoppingData::AddPhysicsVector() for element", "mat040", 
		 FatalException, ed);
    return false;
  }

  dedxMapElements[key] = physicsVector;

  return true;
}

// #########################################################################

G4bool G4IonStoppingData::RemovePhysicsVector(
	G4int atomicNumberIon,         // Atomic number of ion
        const G4String& matIdentifier  // Name or chemical formula of material
				      ) {

  G4IonDEDXKeyMat key = std::make_pair(atomicNumberIon, matIdentifier);

  auto iter = dedxMapMaterials.find(key);

  if(iter == dedxMapMaterials.end()) {
    G4Exception ("G4IonStoppingData::RemovePhysicsVector() for material", 
		 "mat038", FatalException, "Invalid name of the material.");
    return false;
  }

  G4PhysicsVector* physicsVector = (*iter).second;

  // Deleting key of physics vector from material map
  dedxMapMaterials.erase(key);

  // Deleting physics vector
  delete physicsVector;

  return true;
}

// #########################################################################

G4bool G4IonStoppingData::RemovePhysicsVector(
	G4int atomicNumberIon,         // Atomic number of ion
        G4int atomicNumberElem         // Atomic number of elemental material
				      ) {
  G4IonDEDXKeyElem key = std::make_pair(atomicNumberIon, atomicNumberElem);

  auto iter = dedxMapElements.find(key);

  if(iter == dedxMapElements.end()) {
    G4Exception ("G4IonStoppingData::RemovePhysicsVector() for element", 
		 "mat038", FatalException, "Invalid element.");
    return false;
  }

  G4PhysicsVector* physicsVector = (*iter).second;

  // Deleting key of physics vector from material map
  dedxMapElements.erase(key);

  // Deleting physics vector
  delete physicsVector;

  return true;
}

// #########################################################################

G4bool G4IonStoppingData::BuildPhysicsVector(
	G4int atomicNumberIon,          // Atomic number of ion
        const G4String& matname         // Name of material
        					     ) {
  if(IsApplicable(atomicNumberIon, matname))
  {
    return true;
  }

  const char* path = G4FindDataDir("G4LEDATA");
  if(path == nullptr)
  {
    G4Exception("G4IonStoppingData::BuildPhysicsVector()", "mat521",
                FatalException, "G4LEDATA environment variable not set");
    return false;
  }

  std::ostringstream file;
  G4String ww = (fICRU90 && (matname == "G4_WATER" || 
			     matname == "G4_AIR" || 
			     matname == "G4_GRAPHITE")) ? "90" : "73";
 
  file << path << "/" << subDir << ww << "/z" 
       << atomicNumberIon << "_" << matname << ".dat";
  G4String fileName = G4String( file.str().c_str() );

  std::ifstream ifilestream( fileName );

  if(!ifilestream.is_open())
  {
    return false;
  }

  auto* physicsVector = new G4PhysicsFreeVector(true);

  if( !physicsVector -> Retrieve(ifilestream, true) ) {
 
     ifilestream.close();
     return false;
  }

  physicsVector -> ScaleVector( MeV, MeV * cm2 /( 0.001 * g) ); 
  physicsVector -> FillSecondDerivatives();

  // Retrieved vector is added to material store
  if( !AddPhysicsVector(physicsVector, atomicNumberIon, matname) ) {
     delete physicsVector;
     ifilestream.close();
     return false;
  }

  ifilestream.close();
  return true;
}

// #########################################################################

G4bool G4IonStoppingData::BuildPhysicsVector(
	G4int ZIon,          // Atomic number of ion
        G4int ZElem          // Atomic number of elemental material
        					     ) 
{
  if(IsApplicable(ZIon, ZElem))
  {
    return true;
  }

  const char* path = G4FindDataDir("G4LEDATA");
  if(path == nullptr)
  {
    G4Exception("G4IonStoppingData::BuildPhysicsVector()", "mat522",
                FatalException, "G4LEDATA environment variable not set");
    return false;
  }
  std::ostringstream file;
  G4String ww = (fICRU90 && ZIon <= 18 && 
		 (ZElem == 1 || ZElem == 6 || 
		  ZElem == 7 || ZElem == 8)) ? "90" : "73";
 
  file << path << "/" << subDir << ww << "/z" 
       << ZIon << "_" << ZElem << ".dat";
                      
  G4String fileName = G4String( file.str().c_str() );
  std::ifstream ifilestream( fileName );

  if(!ifilestream.is_open())
  {
    return false;
  }
  auto* physicsVector = new G4PhysicsFreeVector(true);

  if( !physicsVector -> Retrieve(ifilestream, true) ) {
     ifilestream.close();
     return false;
  }

  physicsVector -> ScaleVector( MeV, MeV * cm2 /( 0.001 * g) ); 
  physicsVector -> FillSecondDerivatives();

  // Retrieved vector is added to material store
  if( !AddPhysicsVector(physicsVector, ZIon, ZElem) ) {
     delete physicsVector;
     ifilestream.close();
     return false;
  }

  ifilestream.close();
  return true;
}

// #########################################################################

void G4IonStoppingData::ClearTable() {
  auto iterMat     = dedxMapMaterials.begin();
  auto iterMat_end = dedxMapMaterials.end();

  for(;iterMat != iterMat_end; iterMat++) { 

    G4PhysicsVector* vec = iterMat -> second;

    delete vec;
  }

  dedxMapMaterials.clear();

  auto iterElem     = dedxMapElements.begin();
  auto iterElem_end = dedxMapElements.end();

  for(;iterElem != iterElem_end; iterElem++) { 

    G4PhysicsVector* vec = iterElem -> second;

    delete vec;
  }

  dedxMapElements.clear();
}

// #########################################################################

void G4IonStoppingData::DumpMap() {
  auto iterMat     = dedxMapMaterials.begin();
  auto iterMat_end = dedxMapMaterials.end();

  G4cout << std::setw(15) << std::right
         << "Atomic nmb ion"
         << std::setw(25) << std::right
         << "Material name"
         << G4endl;

  for(;iterMat != iterMat_end; iterMat++) {
      G4IonDEDXKeyMat key = iterMat -> first;
      G4PhysicsVector* physicsVector = iterMat -> second; 

      G4int atomicNumberIon = key.first;
      G4String matIdentifier = key.second;

      if(physicsVector != nullptr)
      {
        G4cout << std::setw(15) << std::right << atomicNumberIon
               << std::setw(25) << std::right << matIdentifier << G4endl;
      }
  }

  auto iterElem     = dedxMapElements.begin();
  auto iterElem_end = dedxMapElements.end();

  G4cout << std::setw(15) << std::right
         << "Atomic nmb ion"
         << std::setw(25) << std::right
         << "Atomic nmb material"
         << G4endl;

  for(;iterElem != iterElem_end; iterElem++) { 
      G4IonDEDXKeyElem key = iterElem -> first;
      G4PhysicsVector* physicsVector = iterElem -> second; 

      G4int atomicNumberIon = key.first;
      G4int atomicNumberElem = key.second;

      if(physicsVector != nullptr)
      {
        G4cout << std::setw(15) << std::right << atomicNumberIon
               << std::setw(25) << std::right << atomicNumberElem << G4endl;
      }
  }

}

// #########################################################################

