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
// GEANT4 class header file
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
//
//
// Class description: Class which can read ion stopping power data from
//                    $G4LEDATA/ion_stopping_data
//
// Comments:
//
// =========================================================================== 
//

#ifndef G4IONSTOPPINGDATA_HH
#define G4IONSTOPPINGDATA_HH

#include "globals.hh"
#include "G4VIonDEDXTable.hh"
#include <utility>
#include <vector>
#include <map>


class G4IonStoppingData : public G4VIonDEDXTable {

public:
 G4IonStoppingData(const G4String& dir, G4bool icru);
 ~G4IonStoppingData() override;

 // Function for checking the availability of stopping power tables
 // for a given ion-material couple, where the material consists of
 // a single element only.
 G4bool IsApplicable(
   G4int atomicNumberIon,  // Atomic number of ion
   G4int atomicNumberElem  // Atomic number of elemental material
   ) override;

 // Function for checking the availability of stopping power tables
 // for given ion-material couples.
 G4bool IsApplicable(
   G4int atomicNumberIon,         // Atomic number of ion
   const G4String& matIdentifier  // Name or chemical formula of material
   ) override;

 // Function which invokes the read/build process of physics vectors from
 // files in G4LEDATA
 G4bool BuildPhysicsVector(G4int ionZ, const G4String& matName) override;

 // Function which invokes the read/build process of physics vectors from
 // files in G4LEDATA
 G4bool BuildPhysicsVector(G4int ionZ, G4int matZ) override;

 // Function returning the stopping power vector for given ion-material
 // couples, where the material consists of a single element only.
 G4PhysicsVector* GetPhysicsVector(
   G4int atomicNumberIon,  // Atomic number of ion
   G4int atomicNumberElem  // Atomic number of elemental material
   ) override;

 // Function returning the stopping power vector for given ion-material
 // couples.
 G4PhysicsVector* GetPhysicsVector(
   G4int atomicNumberIon,       // Atomic number of ion
   const G4String& matIdenfier  // Name or chemical formula of material
   ) override;

 // Function returning the stopping power value for given ion-material
 // couples, where the material consists of a single element only, and
 // given energy.
 G4double GetDEDX(G4double kinEnergyPerNucleon,  // Kinetic energy per nucleon
                  G4int atomicNumberIon,         // Atomic number of ion
                  G4int atomicNumberElem  // Atomic number of elemental material
 );

 // Function returning the stopping power value for given ion-material
 // couples and given energy.
 G4double GetDEDX(
   G4double kinEnergyPerNucleon,  // Kinetic energy per nucleon
   G4int atomicNumberIon,         // Atomic number of ion
   const G4String& matIdentifier  // Name or chemical formula of material
 );

 // Function for adding dE/dx vector for an elemental materials. The last
 // argument only applies to elemental materials.
 G4bool AddPhysicsVector(
   G4PhysicsVector* physicsVector,  // Physics vector
   G4int atomicNumberIon,           // Atomic number of ion
   const G4String& matIdentifier    // Name or chemical formula of material
 );

 // Function for adding dE/dx vector for an elemental materials. The last
 // argument only applies to elemental materials.
 G4bool AddPhysicsVector(
   G4PhysicsVector* physicsVector,  // Physics vector
   G4int atomicNumberIon,           // Atomic number of ion
   G4int atomicNumberElem           // Atomic number of elemental material
 );

 // Function for removing dE/dx vector for a compound materials
 G4bool RemovePhysicsVector(
   G4int atomicNumberIon,         // Atomic number of ion
   const G4String& matIdentifier  // Name or chemical formula of material
 );
 // Function for removing dE/dx vector for a compound materials
 G4bool RemovePhysicsVector(
   G4int atomicNumberIon,  // Atomic number of ion
   G4int atomicNumberElem  // Atomic number of elemental material
 );
 // Function deleting all physics vectors and clearing the maps
 void ClearTable();

 // Function printing the ion-material pairs of available vectors to stdout
 void DumpMap();

private:
 G4IonStoppingData(G4IonStoppingData&) = delete;
 const G4IonStoppingData& operator=(const G4IonStoppingData&) = delete;

 // Subdirectory of G4LEDATA
 G4String subDir;

 using G4IonDEDXKeyElem = std::pair<G4int, G4int>;
 using G4IonDEDXKeyMat  = std::pair<G4int, G4String>;

 using G4IonDEDXMapElem = std::map<G4IonDEDXKeyElem, G4PhysicsVector*>;
 using G4IonDEDXMapMat  = std::map<G4IonDEDXKeyMat, G4PhysicsVector*>;

 G4IonDEDXMapElem dedxMapElements;
 G4IonDEDXMapMat dedxMapMaterials;

 G4bool fICRU90;
};

#endif // G4IONSTOPPINGDATA_HH
