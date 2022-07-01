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
// Class:                G4ExtDEDXTable
//
// Base class:           G4VIonDEDXTable 
// 
// Author:               Anton Lechner (Anton.Lechner@cern.ch)
//
// First implementation: 29. 02. 2009
//
// Modifications:
// 03.11.2009 A. Lechner:  Added new methods BuildPhysicsVector according
//            to interface changes in base class G4VIonDEDXTable.
// 25.10.2010 V.Ivanchenko fixed bug in usage of iterators reported by the 
//            Coverity tool
// 01.11.2010 V.Ivanchenko fixed remaining bugs reported by Coverity 
//
//
// Class description:
//    Utility class for users to add their own electronic stopping powers
//    for ions. This class is dedicated for use with G4IonParametrisedLossModel
//    of the low-energy electromagnetic package.
//
// Comments:
//
// =========================================================================== 
//

#include "G4ExtDEDXTable.hh" 
#include "G4PhysicsVector.hh"
#include "G4PhysicsVectorType.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
#include <fstream>
#include <sstream>
#include <iomanip>

// #########################################################################

G4ExtDEDXTable::~G4ExtDEDXTable() {

  ClearTable();
}

// #########################################################################

G4bool G4ExtDEDXTable::BuildPhysicsVector(G4int ionZ, G4int matZ) {

  return IsApplicable( ionZ, matZ );
}


// #########################################################################

G4bool G4ExtDEDXTable::BuildPhysicsVector(G4int ionZ, 
                                          const G4String& matName) {

  return IsApplicable( ionZ, matName );
}

// #########################################################################

G4bool G4ExtDEDXTable::IsApplicable(
         G4int atomicNumberIon,  // Atomic number of ion
         G4int atomicNumberElem  // Atomic number of elemental material
				    )
{
  G4IonDEDXKeyElem key = std::make_pair(atomicNumberIon, atomicNumberElem);

  auto iter = dedxMapElements.find(key);

  return iter != dedxMapElements.end();
}

// #########################################################################

G4bool G4ExtDEDXTable::IsApplicable(
         G4int atomicNumberIon,         // Atomic number of ion
         const G4String& matIdentifier  // Name or chemical formula of material
				    )
{
  G4IonDEDXKeyMat key = std::make_pair(atomicNumberIon, matIdentifier);

  auto iter = dedxMapMaterials.find(key);

  return iter != dedxMapMaterials.end();
}

// #########################################################################

G4PhysicsVector* G4ExtDEDXTable::GetPhysicsVector(
         G4int atomicNumberIon,        // Atomic number of ion
         G4int atomicNumberElem        // Atomic number of elemental material
				    )
{
  G4IonDEDXKeyElem key = std::make_pair(atomicNumberIon, atomicNumberElem);

  auto iter = dedxMapElements.find(key);

  return (iter != dedxMapElements.end()) ? iter->second : nullptr; 
}

// #########################################################################

G4PhysicsVector*  G4ExtDEDXTable::GetPhysicsVector(
         G4int atomicNumberIon,        // Atomic number of ion
         const G4String& matIdentifier // Name or chemical formula of material
				    )
{
  G4IonDEDXKeyMat key = std::make_pair(atomicNumberIon, matIdentifier);

  auto iter = dedxMapMaterials.find(key);

  return (iter != dedxMapMaterials.end()) ? iter->second : nullptr; 
}

// #########################################################################

G4double G4ExtDEDXTable::GetDEDX(
         G4double kinEnergyPerNucleon, // Kinetic energy per nucleon
         G4int atomicNumberIon,        // Atomic number of ion
         G4int atomicNumberElem        // Atomic number of elemental material
				  )
{
  G4IonDEDXKeyElem key = std::make_pair(atomicNumberIon, atomicNumberElem);

  auto iter = dedxMapElements.find(key);

  return ( iter != dedxMapElements.end() ) ?
    (iter->second)->Value( kinEnergyPerNucleon) : 0.0;
}

// #########################################################################

G4double G4ExtDEDXTable::GetDEDX(
         G4double kinEnergyPerNucleon, // Kinetic energy per nucleon
         G4int atomicNumberIon,        // Atomic number of ion
         const G4String& matIdentifier // Name or chemical formula of material
				  )
{
  G4IonDEDXKeyMat key = std::make_pair(atomicNumberIon, matIdentifier);

  auto iter = dedxMapMaterials.find(key);

  return (iter != dedxMapMaterials.end()) ?
    (iter->second)->Value( kinEnergyPerNucleon) : 0.0;
}

// #########################################################################

G4bool G4ExtDEDXTable::AddPhysicsVector(
        G4PhysicsVector* physicsVector, // Physics vector
	G4int atomicNumberIon,          // Atomic number of ion
        const G4String& matIdentifier,  // Name of elemental material
        G4int atomicNumberElem          // Atomic number of elemental material
				      ) {

  if(physicsVector == nullptr) {
    G4Exception ("G4ExtDEDXTable::AddPhysicsVector() for material", 
		 "mat037", FatalException, 
		 "Pointer to vector is null-pointer.");
    return false;
  }

  if(matIdentifier.empty()) {
    G4Exception ("G4ExtDEDXTable::AddPhysicsVector() for material", 
                 "mat038", FatalException, "Invalid name of the material.");
     return false;
  }

  if(atomicNumberIon <= 2) {
    G4Exception ("G4ExtDEDXTable::AddPhysicsVector() for material", 
                 "mat039", FatalException, "Illegal atomic number.");
    return false;
  }

  if(atomicNumberElem > 0) {

     G4IonDEDXKeyElem key = std::make_pair(atomicNumberIon, atomicNumberElem);

     if(dedxMapElements.count(key) == 1) {
       G4Exception ("G4ExtDEDXTable::AddPhysicsVector() for material", 
		    "mat037", FatalException, 
		    "Vector already exist, remove it before replacing.");
       return false;
     }

     dedxMapElements[key] = physicsVector;
  }

  G4IonDEDXKeyMat mkey = std::make_pair(atomicNumberIon, matIdentifier);

  if(dedxMapMaterials.count(mkey) == 1) {
    G4Exception ("G4ExtDEDXTable::AddPhysicsVector() for material", 
		 "mat037", FatalException, 
		 "Vector already exist, remove it before replacing.");
    return false;
  }

  dedxMapMaterials[mkey] = physicsVector;

  return true;
}

// #########################################################################

G4bool G4ExtDEDXTable::RemovePhysicsVector(
	G4int atomicNumberIon,         // Atomic number of ion
        const G4String& matIdentifier  // Name or chemical formula of material
				      ) {

  G4PhysicsVector* physicsVector = nullptr;

  // Deleting key of physics vector from material map
  G4IonDEDXKeyMat key = std::make_pair(atomicNumberIon, matIdentifier);

  auto iter = dedxMapMaterials.find(key);

  if(iter == dedxMapMaterials.end()) {
    G4Exception ("G4ExtDEDXTable::RemovePhysicsVector() for material", 
		 "mat037", FatalException, 
		 "Pointer to vector is null-pointer.");
    return false;
  }

  physicsVector = (*iter).second;
  dedxMapMaterials.erase(key);

  // Deleting key of physics vector from elemental material map (if it exists)
  G4IonDEDXMapElem::iterator it;
  
  for(it=dedxMapElements.begin(); it!=dedxMapElements.end(); ++it) {

     if( (*it).second == physicsVector ) {
        dedxMapElements.erase(it);
        break;
     }
  }

  // Deleting physics vector
  delete physicsVector;

  return true;
}

// #########################################################################

G4bool G4ExtDEDXTable::StorePhysicsTable(
         const G4String& fileName // File name
				    ) {
  G4bool success = true;

  std::ofstream ofilestream;

  ofilestream.open( fileName, std::ios::out );

  if( !ofilestream ) {
    G4ExceptionDescription ed;
    ed << "Cannot open file " << fileName; 
    G4Exception ("G4IonStoppingData::StorePhysicsTable()", 
                 "mat030", FatalException, ed);
    success = false;
  }   
  else {

     size_t nmbMatTables = dedxMapMaterials.size();

     ofilestream << nmbMatTables << G4endl << G4endl;

     auto iterMat     = dedxMapMaterials.begin();
     auto iterMat_end = dedxMapMaterials.end();

     for(;iterMat != iterMat_end; iterMat++) {
         G4IonDEDXKeyMat key = iterMat -> first;
         G4PhysicsVector* physicsVector = iterMat -> second; 

         G4int atomicNumberIon = key.first;
         G4String matIdentifier = key.second;

         G4int atomicNumberElem = FindAtomicNumberElement(physicsVector);

         if(physicsVector != nullptr) {
  	    ofilestream << atomicNumberIon << "  " << matIdentifier;

        if(atomicNumberElem > 0)
        {
          ofilestream << "  " << atomicNumberElem;
        }

            ofilestream << "  # <Atomic number ion>  <Material name>  ";

            if(atomicNumberElem > 0)
            {
              ofilestream << "<Atomic number element>";
            }

            ofilestream << G4endl << physicsVector -> GetType() << G4endl;

            physicsVector -> Store(ofilestream, true);

            ofilestream << G4endl;
         } else {
	   G4Exception ("G4IonStoppingData::StorePhysicsTable()", 
			"mat030", FatalException,"Cannot store vector.");
         }
     }
  }

  ofilestream.close();

  return success; 
}

// #########################################################################

G4bool G4ExtDEDXTable::RetrievePhysicsTable(const G4String& fileName) 
{ 
  std::ifstream ifilestream;
  ifilestream.open( fileName, std::ios::in|std::ios::binary );
  if( ! ifilestream ) {
    G4ExceptionDescription ed;
    ed << "Cannot open file " << fileName; 
    G4Exception ("G4IonStoppingData::RetrievePhysicsTable()", 
                 "mat030", FatalException, ed);
    return false;
  }   

  //std::string::size_type nmbVectors;
  G4int nmbVectors = 0;
  ifilestream >> nmbVectors;
  if( ifilestream.fail() || nmbVectors <= 0) { 
    G4cout << "G4ExtDEDXTable::RetrievePhysicsTable() " 
	   << " File content of " << fileName << " ill-formated."
	   << " Nvectors= " << nmbVectors
	   << G4endl;
    ifilestream.close(); 
    return false; 
  }

  for(G4int i = 0; i<nmbVectors; ++i) {

    G4String line = "";
    // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    while( line.empty() ) {

      getline( ifilestream, line );
      if( ifilestream.fail() ) { 
	G4cout << "G4ExtDEDXTable::RetrievePhysicsTable() " 
	       << " File content of " << fileName << " ill-formated." 
	       << G4endl;
	ifilestream.close(); 
	return false; 
      }

      std::string::size_type pos = line.find_first_of("#");
      if(pos != std::string::npos && pos > 0) {
	line = line.substr(0, pos);
      }
    }

    std::istringstream headerstream( line );     

    std::string::size_type atomicNumberIon;
    headerstream >> atomicNumberIon;

    G4String materialName;
    headerstream >> materialName;

    if( headerstream.fail() || std::string::npos == atomicNumberIon) {
      G4cout << "G4ExtDEDXTable::RetrievePhysicsTable() " 
	     << " File content of " << fileName << " ill-formated "
	     << " (vector header)." 
	     << G4endl;
      ifilestream.close();
      return false;
    } 

    std::string::size_type atomicNumberMat;
    headerstream >> atomicNumberMat;

    if( headerstream.eof() || std::string::npos == atomicNumberMat) { 
      atomicNumberMat = 0; 
    }

    G4int vectorType;
    ifilestream >> vectorType;
      
    G4PhysicsVector* physicsVector = CreatePhysicsVector(vectorType);

    if(physicsVector == nullptr) {
      G4cout << "G4ExtDEDXTable::RetrievePhysicsTable  "
	     << " illegal physics Vector type " << vectorType
	     << " in  " << fileName 
	     << G4endl;
      ifilestream.close();
      return false;
    }

    if( !physicsVector -> Retrieve(ifilestream, true) ) {
      G4cout << "G4ExtDEDXTable::RetrievePhysicsTable() " 
	     << " File content of " << fileName << " ill-formated." 
	     << G4endl;
      ifilestream.close();
      return false;
    } 
    physicsVector -> FillSecondDerivatives();

    // Retrieved vector is added to material store
    if( !AddPhysicsVector(physicsVector, (G4int)atomicNumberIon, 
			  materialName, (G4int)atomicNumberMat) ) {

      delete physicsVector;
      ifilestream.close();
      return false;
    }
  }

  ifilestream.close();

  return true;
}

// #########################################################################

G4PhysicsVector* G4ExtDEDXTable::CreatePhysicsVector(G4int vectorType) {

  G4PhysicsVector* physicsVector = nullptr;

  switch (vectorType) {

  case T_G4PhysicsLinearVector: 
    physicsVector = new G4PhysicsLinearVector(true);
    break;

  case T_G4PhysicsLogVector: 
    physicsVector = new G4PhysicsLogVector(true);
    break;

  case T_G4PhysicsFreeVector: 
    physicsVector = new G4PhysicsFreeVector(true);
    break;
  
  default:
    break;
  }
  return physicsVector;
}

// #########################################################################

G4int G4ExtDEDXTable::FindAtomicNumberElement(
                                   G4PhysicsVector* physicsVector
                                                   ) {

  G4int atomicNumber = 0;

  auto iter     = dedxMapElements.begin();
  auto iter_end = dedxMapElements.end();

  for(;iter != iter_end; ++iter) {

     if( (*iter).second == physicsVector ) {

        G4IonDEDXKeyElem key = (*iter).first;
        atomicNumber = key.second;
     }
  }

  return atomicNumber;
}

// #########################################################################

void G4ExtDEDXTable::ClearTable() {
  auto iterMat     = dedxMapMaterials.begin();
  auto iterMat_end = dedxMapMaterials.end();

  for(;iterMat != iterMat_end; ++iterMat) { 

    G4PhysicsVector* vec = iterMat -> second;

    delete vec;
  }

  dedxMapElements.clear();
  dedxMapMaterials.clear();
}

// #########################################################################

void G4ExtDEDXTable::DumpMap() {
  auto iterMat     = dedxMapMaterials.begin();
  auto iterMat_end = dedxMapMaterials.end();

  G4cout << std::setw(15) << std::right
         << "Atomic nmb ion"
         << std::setw(25) << std::right
         << "Material name"
         << std::setw(25) << std::right
         << "Atomic nmb material"
         << G4endl;

  for(;iterMat != iterMat_end; ++iterMat) {
      G4IonDEDXKeyMat key = iterMat -> first;
      G4PhysicsVector* physicsVector = iterMat -> second; 

      G4int atomicNumberIon = key.first;
      G4String matIdentifier = key.second;

      G4int atomicNumberElem = FindAtomicNumberElement(physicsVector);

      if(physicsVector != nullptr)
      {
        G4cout << std::setw(15) << std::right << atomicNumberIon
               << std::setw(25) << std::right << matIdentifier << std::setw(25)
               << std::right;

        if(atomicNumberElem > 0)
        {
          G4cout << atomicNumberElem;
        }
        else
        {
          G4cout << "N/A";
        }

        G4cout << G4endl;
      }
  }

}

// #########################################################################

