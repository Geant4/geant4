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
// $Id: G4IonStoppingData.cc,v 1.3 2010-10-25 08:41:39 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//
//
// Class description: Class which can read ion stopping power data from
//                    $G4LEDATA/ion_stopping_data
//
// Comments:
//
// =========================================================================== 
//

#include "G4IonStoppingData.hh" 
#include "G4PhysicsVector.hh"
#include "G4LPhysicsFreeVector.hh"
#include <fstream>
#include <sstream>
#include <iomanip>


// #########################################################################

G4IonStoppingData::G4IonStoppingData(const G4String& leDirectory) :
  subDir( leDirectory ) {

}

// #########################################################################

G4IonStoppingData::~G4IonStoppingData() {

  ClearTable();
}

// #########################################################################

G4bool G4IonStoppingData::IsApplicable(
         G4int atomicNumberIon,  // Atomic number of ion
         G4int atomicNumberElem  // Atomic number of elemental material
				    ) {
  G4bool isApplicable = true; 
  G4IonDEDXKeyElem key = std::make_pair(atomicNumberIon, atomicNumberElem);

  G4IonDEDXMapElem::iterator iter = dedxMapElements.find(key);

  if(iter == dedxMapElements.end()) isApplicable = false; 

  return isApplicable; 
}

// #########################################################################

G4bool G4IonStoppingData::IsApplicable(
         G4int atomicNumberIon,         // Atomic number of ion
         const G4String& matIdentifier  // Name or chemical formula of material
				    ) {
  G4bool isApplicable = true; 
  G4IonDEDXKeyMat key = std::make_pair(atomicNumberIon, matIdentifier);

  G4IonDEDXMapMat::iterator iter = dedxMapMaterials.find(key);

  if(iter == dedxMapMaterials.end()) isApplicable = false; 

  return isApplicable; 
}

// #########################################################################

G4PhysicsVector* G4IonStoppingData::GetPhysicsVector(
         G4int atomicNumberIon,        // Atomic number of ion
         G4int atomicNumberElem        // Atomic number of elemental material
				    ) {

  G4PhysicsVector* physVector = 0;

  G4IonDEDXKeyElem key = std::make_pair(atomicNumberIon, atomicNumberElem);

  G4IonDEDXMapElem::iterator iter = dedxMapElements.find(key);

  if(iter != dedxMapElements.end()) physVector = iter -> second; 

  return physVector; 
}

// #########################################################################

G4PhysicsVector*  G4IonStoppingData::GetPhysicsVector(
         G4int atomicNumberIon,        // Atomic number of ion
         const G4String& matIdentifier // Name or chemical formula of material
				    ) {

  G4PhysicsVector* physVector = 0;

  G4IonDEDXKeyMat key = std::make_pair(atomicNumberIon, matIdentifier);

  G4IonDEDXMapMat::iterator iter = dedxMapMaterials.find(key);

  if(iter != dedxMapMaterials.end()) physVector = iter -> second; 

  return physVector; 
}

// #########################################################################

G4double G4IonStoppingData::GetDEDX(
         G4double kinEnergyPerNucleon, // Kinetic energy per nucleon
         G4int atomicNumberIon,        // Atomic number of ion
         G4int atomicNumberElem        // Atomic number of elemental material
				  ) {
  G4double dedx = 0;

  G4IonDEDXKeyElem key = std::make_pair(atomicNumberIon, atomicNumberElem);

  G4IonDEDXMapElem::iterator iter = dedxMapElements.find(key);

  if( iter != dedxMapElements.end() ) {
     G4PhysicsVector* physicsVector = iter -> second; 

     G4bool b;
     dedx = physicsVector -> GetValue( kinEnergyPerNucleon, b );   
  }

  return dedx; 
}

// #########################################################################

G4double G4IonStoppingData::GetDEDX(
         G4double kinEnergyPerNucleon, // Kinetic energy per nucleon
         G4int atomicNumberIon,        // Atomic number of ion
         const G4String& matIdentifier // Name or chemical formula of material
				  ) {
  G4double dedx = 0;

  G4IonDEDXKeyMat key = std::make_pair(atomicNumberIon, matIdentifier);

  G4IonDEDXMapMat::iterator iter = dedxMapMaterials.find(key);

  if(iter != dedxMapMaterials.end()) {
     G4PhysicsVector* physicsVector = iter -> second; 

     G4bool b;
     dedx = physicsVector -> GetValue( kinEnergyPerNucleon, b );   
  }

  return dedx; 
}

// #########################################################################

G4bool G4IonStoppingData::AddPhysicsVector(
        G4PhysicsVector* physicsVector, // Physics vector
	G4int atomicNumberIon,          // Atomic number of ion
        const G4String& matIdentifier   // Name of elemental material
				      ) {

  if(physicsVector == 0) {

#ifdef G4VERBOSE
     G4cerr << "G4IonStoppingData::AddPhysicsVector() Error: Pointer to vector"
            << " is null-pointer."
            << G4endl;
#endif

     return false;
  }

  if(matIdentifier.empty()) {

#ifdef G4VERBOSE
     G4cerr << "G4IonStoppingData::AddPhysicsVector() Error: "
            << "Cannot add physics vector. Invalid name."
            << G4endl;
#endif

     return false;
  }

  if(atomicNumberIon <= 0) {

#ifdef G4VERBOSE
     G4cerr << "G4IonStoppingData::AddPhysicsVector() Error: "
            << "Cannot add physics vector. Illegal atomic number."
            << G4endl;
#endif

     return false;
  }

  G4IonDEDXKeyMat mkey = std::make_pair(atomicNumberIon, matIdentifier);

  if(dedxMapMaterials.count(mkey) == 1) {

#ifdef G4VERBOSE
     G4cerr << "G4IonStoppingData::AddPhysicsVector() Error: "
            << "Vector with Z1 = " << atomicNumberIon << ", mat = " 
            << matIdentifier
            << "already exists. Remove first before replacing."
            << G4endl;
#endif

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
				      ) {

  if(physicsVector == 0) {

#ifdef G4VERBOSE
     G4cerr << "G4IonStoppingData::AddPhysicsVector() Error: "
            << "Pointer to vector is null-pointer."
            << G4endl;
#endif

     return false;
  }

  if(atomicNumberIon <= 0) {

#ifdef G4VERBOSE
     G4cerr << "G4IonStoppingData::AddPhysicsVector() Error: "
            << "Cannot add physics vector. Illegal atomic number."
            << G4endl;
#endif

     return false;
  }

  if(atomicNumberElem <= 0) {

#ifdef G4VERBOSE
        G4cerr << "G4IonStoppingData::AddPhysicsVector() Error: "
               << "Atomic number of element < 0."
               << G4endl;
#endif
        return false;
  }

  G4IonDEDXKeyElem key = std::make_pair(atomicNumberIon, atomicNumberElem);

  if(dedxMapElements.count(key) == 1) {

#ifdef G4VERBOSE
     G4cerr << "G4IonStoppingData::AddPhysicsVector() Error: "
            << "Vector with Z1 = " << atomicNumberIon << ", Z2 = " 
            << atomicNumberElem
            << " already exists. Remove first before replacing."
            << G4endl;
#endif
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

  G4IonDEDXMapMat::iterator iter = dedxMapMaterials.find(key);

  if(iter == dedxMapMaterials.end()) {

#ifdef G4VERBOSE
     G4cerr << "G4IonStoppingData::RemovePhysicsVector() Warning: "
            << "Cannot remove physics vector. Vector not found."
            << G4endl;
#endif

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

  G4IonDEDXMapElem::iterator iter = dedxMapElements.find(key);

  if(iter == dedxMapElements.end()) {

#ifdef G4VERBOSE
     G4cerr << "G4IonStoppingData::RemovePhysicsVector() Warning: "
            << "Cannot remove physics vector. Vector not found."
            << G4endl;
#endif

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
        const G4String& matIdentifier   // Name of material
        					     ) {

  if( IsApplicable(atomicNumberIon, matIdentifier) ) return true;

  char* path = getenv("G4LEDATA");
  if ( !path ) {
    G4Exception("G4IonStoppingData::BuildPhysicsVector: G4LEDATA environment variable not set");
    return false;
  }
  
  std::ostringstream file;
 
  file << path << "/" << subDir << "/z" 
       << atomicNumberIon << "_" << matIdentifier 
       << ".dat";
                      
  G4String fileName = G4String( file.str().c_str() );

  std::ifstream ifilestream( fileName );

  if ( !ifilestream.is_open() ) return false;
  
  G4LPhysicsFreeVector* physicsVector = new G4LPhysicsFreeVector(); 

  if( !physicsVector -> Retrieve(ifilestream, true) ) {
 
     ifilestream.close();
     return false;
  }

  physicsVector -> ScaleVector( MeV, MeV * cm2 /( 0.001 * g) ); 
  physicsVector -> SetSpline( true );
  physicsVector -> FillSecondDerivatives();

  // Retrieved vector is added to material store
  if( !AddPhysicsVector(physicsVector, atomicNumberIon, matIdentifier) ) {
     delete physicsVector;
     ifilestream.close();
     return false;
  }

  ifilestream.close();
  return true;
}

// #########################################################################

G4bool G4IonStoppingData::BuildPhysicsVector(
	G4int atomicNumberIon,          // Atomic number of ion
        G4int atomicNumberElem          // Atomic number of elemental material
        					     ) {

  if( IsApplicable(atomicNumberIon, atomicNumberElem) ) return true;

  char* path = getenv("G4LEDATA");
  if ( !path ) {
    G4Exception("G4IonStoppingData::BuildPhysicsVector: G4LEDATA environment variable not set");
    return false;
  }
  std::ostringstream file;
 
  file << path << "/" << subDir << "/z" 
       << atomicNumberIon << "_" << atomicNumberElem
       << ".dat";
                      
  G4String fileName = G4String( file.str().c_str() );

  std::ifstream ifilestream( fileName );

  if ( !ifilestream.is_open() ) return false;
  
  G4LPhysicsFreeVector* physicsVector = new G4LPhysicsFreeVector(); 

  if( !physicsVector -> Retrieve(ifilestream, true) ) {
 
     ifilestream.close();
     return false;
  }

  physicsVector -> ScaleVector( MeV, MeV * cm2 /( 0.001 * g) ); 
  physicsVector -> SetSpline( true );
  physicsVector -> FillSecondDerivatives();

  // Retrieved vector is added to material store
  if( !AddPhysicsVector(physicsVector, atomicNumberIon, atomicNumberElem) ) {
     delete physicsVector;
     ifilestream.close();
     return false;
  }

  ifilestream.close();
  return true;
}

// #########################################################################

void G4IonStoppingData::ClearTable() {

  G4IonDEDXMapMat::iterator iterMat = dedxMapMaterials.begin();
  G4IonDEDXMapMat::iterator iterMat_end = dedxMapMaterials.end();

  for(;iterMat != iterMat_end; iterMat++) { 

    G4PhysicsVector* vec = iterMat -> second;

    if(vec != 0) delete vec;
  }

  dedxMapMaterials.clear();

  G4IonDEDXMapElem::iterator iterElem = dedxMapElements.begin();
  G4IonDEDXMapElem::iterator iterElem_end = dedxMapElements.end();

  for(;iterElem != iterElem_end; iterElem++) { 

    G4PhysicsVector* vec = iterElem -> second;

    if(vec != 0) delete vec;
  }

  dedxMapElements.clear();
}

// #########################################################################

void G4IonStoppingData::DumpMap() {

  G4IonDEDXMapMat::iterator iterMat = dedxMapMaterials.begin();
  G4IonDEDXMapMat::iterator iterMat_end = dedxMapMaterials.end();

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

      if(physicsVector != 0) {
         G4cout << std::setw(15) << std::right
                << atomicNumberIon
                << std::setw(25) << std::right
                << matIdentifier
                << G4endl;
      }
  }

  G4IonDEDXMapElem::iterator iterElem = dedxMapElements.begin();
  G4IonDEDXMapElem::iterator iterElem_end = dedxMapElements.end();

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

      if(physicsVector != 0) {
         G4cout << std::setw(15) << std::right
                << atomicNumberIon
                << std::setw(25) << std::right
                << atomicNumberElem
                << G4endl;
      }
  }

}

// #########################################################################

