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
//
// ===========================================================================
// GEANT4 class header file
//
// Class:                G4VIonDEDXTable 
// 
// Author:               Anton Lechner (Anton.Lechner@cern.ch)
//
// First implementation: 15. 02. 2009
//
// Modifications:        
// 02. 11. 2009 A. Lechner: Added BuildPhysicsVector function which should
//              invoke the build process of physics vectors.
//
// Class description:
//    Base class for electronic stopping power tables for ions.
//
// Comments:
//
// =========================================================================== 
//

#ifndef G4VIONDEDXTABLE_HH
#define G4VIONDEDXTABLE_HH

#include "globals.hh"

class G4PhysicsVector;


class G4VIonDEDXTable {

 public:
   explicit G4VIonDEDXTable() = default;
   virtual ~G4VIonDEDXTable() = default;

   // Function for building a stopping power table for a material consisting of
   // a single element only.
   virtual G4bool BuildPhysicsVector(
        G4int,  // Atomic number of ion
        G4int   // Atomic number of elemental material
                               ) = 0;

   // Function building stopping power table for a given ion-material pair.
   virtual G4bool BuildPhysicsVector(
        G4int,           // Atomic number of ion
        const G4String&  // Name or chemical formula of material
                               ) = 0;

   // Function for checking the availability of stopping power tables
   // for a given ion-material pair, where the material consists of
   // a single element only.
   virtual G4bool IsApplicable(
        G4int,  // Atomic number of ion
        G4int   // Atomic number of elemental material
                               ) = 0;

   // Function for checking the availability of stopping power tables
   // for a given ion-material pair.
   virtual G4bool IsApplicable(
        G4int,           // Atomic number of ion
        const G4String&  // Name or chemical formula of material
                               ) = 0;

   // Function returning the stopping power vector for a given ion-material
   // pair, where the material consists of a single element only (the 
   // physics vector returned by this function must contain the kinetic 
   // energy per nucleon vs the mass stopping power of the elemental material). 
   virtual G4PhysicsVector* GetPhysicsVector(
	G4int,    // Atomic number of ion
        G4int     // Atomic number of elemental material
					     ) = 0;

   // Function returning the stopping power vector for given ion-material
   // pair (the physics vector returned by this function must contain the 
   // kinetic energy per nucleon vs the mass stopping power of the material). 
   virtual G4PhysicsVector* GetPhysicsVector(
	G4int,           // Atomic number of ion
        const G4String&  // Name or chemical formula of material
					     ) = 0;

   G4VIonDEDXTable(G4VIonDEDXTable&) = delete;
   const G4VIonDEDXTable & operator=(const G4VIonDEDXTable&) = delete;

};

#endif // G4VIONDEDXTABLE_HH
