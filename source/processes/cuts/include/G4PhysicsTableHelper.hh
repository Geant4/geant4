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
//
// $Id: G4PhysicsTableHelper.hh,v 1.2 2004/11/07 01:38:17 kurasige Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
   static G4int verboseLevel;
   // controle flag for output message
};

inline
 void  G4PhysicsTableHelper::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline
 G4int G4PhysicsTableHelper::GetVerboseLevel()
{
  return verboseLevel;
}

#endif
