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
// $Id: G4PhysicsTableHelper.hh 70369 2013-05-29 14:59:24Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// Class Description
//  G4PhysicsTableHelper is a static utility class 
//  for helping proceeses to build their physics table
//
// ------------------------------------------------------------
//   First Implementation          20 Aug. 2004   H.Kurashige
//
// ------------------------------------------------------------

#ifndef G4PhysicsTableHelper_h 
#define G4PhysicsTableHelper_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <vector>
#include "G4PhysicsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Region.hh"

class G4PhysicsTableHelper
{ 
  protected:
    G4PhysicsTableHelper();
    ~G4PhysicsTableHelper();
    G4PhysicsTableHelper(const G4PhysicsTableHelper& right);
    G4PhysicsTableHelper& operator=(const G4PhysicsTableHelper&);
 
  public: // with description
    static G4PhysicsTable* PreparePhysicsTable(G4PhysicsTable* physTable);
    // Prepare the given physics table before building the physics table 
    // resize the given physics table to match with the current
    // production cut table.   

    static G4bool RetrievePhysicsTable(G4PhysicsTable* physTable,
				       const G4String& fileName,
				       G4bool ascii              );
    // Retrieve physics table from the given file and 
    // fill the given physics table with retrievd physics vectors

    static void SetPhysicsVector(G4PhysicsTable* physTable,
				 size_t idx,
				 G4PhysicsVector* vec);
    // Set a physics vector at given position 

  public: // with description
   static void  SetVerboseLevel(G4int value);
   static G4int GetVerboseLevel();
   // set/get controle flag for output message
   //  0: Silent
   //  1: Warning message
   //  2: More


 protected:
   static G4ThreadLocal G4int verboseLevel;
   // controle flag for output message
};

#endif
