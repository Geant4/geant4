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

#ifndef G4EXTDEDXTABLE_HH
#define G4EXTDEDXTABLE_HH

#include "globals.hh"
#include "G4VIonDEDXTable.hh"
#include <utility>
#include <vector>
#include <map>


class G4ExtDEDXTable : public G4VIonDEDXTable {

 public:
   explicit G4ExtDEDXTable();
   virtual ~G4ExtDEDXTable();

   virtual G4bool BuildPhysicsVector(G4int ionZ, 
				     const G4String& matName);

   virtual G4bool BuildPhysicsVector(G4int ionZ, 
				     G4int matZ);

   // Function for checking the availability of stopping power tables
   // for a given ion-material couple, where the material consists of
   // a single element only.
   virtual G4bool IsApplicable(
        G4int atomicNumberIon,          // Atomic number of ion
        G4int atomicNumberElem          // Atomic number of elemental material
                       );

   // Function for checking the availability of stopping power tables
   // for given ion-material couples.
   virtual G4bool IsApplicable(
        G4int atomicNumberIon,          // Atomic number of ion
        const G4String& matIdentifier   // Name or chemical formula of material
                       );

   // Function returning the stopping power vector for given ion-material
   // couples, where the material consists of a single element only.
   virtual G4PhysicsVector* GetPhysicsVector(
	G4int atomicNumberIon,          // Atomic number of ion
        G4int atomicNumberElem          // Atomic number of elemental material
				     );

   // Function returning the stopping power vector for given ion-material
   // couples.
   virtual G4PhysicsVector* GetPhysicsVector(
	G4int atomicNumberIon,          // Atomic number of ion
        const G4String& matIdenfier     // Name or chemical formula of material
				     );

   // Function returning the stopping power value for given ion-material
   // couples, where the material consists of a single element only, and
   // given energy.
   G4double GetDEDX(
        G4double kinEnergyPerNucleon,   // Kinetic energy per nucleon
        G4int atomicNumberIon,          // Atomic number of ion
        G4int atomicNumberElem          // Atomic number of elemental material
				     );

   // Function returning the stopping power value for given ion-material
   // couples and given energy.
   G4double GetDEDX(
        G4double kinEnergyPerNucleon,   // Kinetic energy per nucleon
	G4int atomicNumberIon,          // Atomic number of ion
        const G4String& matIdenfier     // Name or chemical formula of material
				     );

   // Function for adding dE/dx vector for an elemental materials. The last
   // argument only applies to elemental materials.
   G4bool AddPhysicsVector(
        G4PhysicsVector* physicsVector, // Physics vector
	G4int atomicNumberIon,          // Atomic number of ion
        const G4String& matIdenfier,    // Name or chemical formula of material
        G4int atomicNumberElem = 0      // Atomic number of elemental material
			 );

   // Function for removing dE/dx vector for a compound materials
   G4bool RemovePhysicsVector(
	G4int atomicNumberIon,          // Atomic number of ion
        const G4String& matIdentifier   // Name or chemical formula of material
			    );

   // Function writing all stopping power vectors to file
   G4bool StorePhysicsTable(
        const G4String& fileName        // File name
			     );

   // Function retrieving all stopping power vectors from file
   G4bool RetrievePhysicsTable(
        const G4String& fileName       // File name
			       );

   // Function deleting all physics vectors and clearing the maps
   void ClearTable();

   // Function printing the ion-material pairs of available vectors to stdout
   void DumpMap();

 private:

   G4ExtDEDXTable(G4ExtDEDXTable&) = delete;
   const G4ExtDEDXTable & operator=(const G4ExtDEDXTable&) = delete;

   G4PhysicsVector* CreatePhysicsVector(G4int vectorType); 

   G4int FindAtomicNumberElement(G4PhysicsVector* physicsVector);

   typedef std::pair<G4int, G4int> G4IonDEDXKeyElem;
   typedef std::pair<G4int, G4String> G4IonDEDXKeyMat;
 
   typedef std::map<G4IonDEDXKeyElem, G4PhysicsVector*> G4IonDEDXMapElem;
   typedef std::map<G4IonDEDXKeyMat, G4PhysicsVector*> G4IonDEDXMapMat;

   G4IonDEDXMapElem dedxMapElements; 
   G4IonDEDXMapMat dedxMapMaterials;
};

#endif // G4EXTDEDXTABLE_HH
