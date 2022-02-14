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
// G4PhysicsTableHelper
//
// Class description:
//
// G4PhysicsTableHelper is a static utility class for helping processes
// to build their physics table.

// Author: H.Kurashige, 20 August 2004 - First implementation
// --------------------------------------------------------------------
#ifndef G4PhysicsTableHelper_hh 
#define G4PhysicsTableHelper_hh 1

#include "globals.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsTable.hh"

class G4PhysicsTableHelper
{ 
  public:

    static G4PhysicsTable* PreparePhysicsTable(G4PhysicsTable* physTable);
      // Prepare the given physics table. Before building the table 
      // resize the given physics table to match with the current
      // production cut table

    static G4bool RetrievePhysicsTable(G4PhysicsTable* physTable,
                                       const G4String& fileName,
                                       G4bool ascii, G4bool spline);
      // Retrieve the physics table from the given file and 
      // fill the given physics table with retrieved physics vectors

    static void SetPhysicsVector(G4PhysicsTable* physTable,
                                 std::size_t idx,
                                 G4PhysicsVector* vec);
      // Set a physics vector at given position 

    static void  SetVerboseLevel(G4int value);
    static G4int GetVerboseLevel();
      // Set/get control flag for output message
      //  0: Silent
      //  1: Warning message
      //  2: More

    G4PhysicsTableHelper(const G4PhysicsTableHelper&) = delete;
    G4PhysicsTableHelper& operator=(const G4PhysicsTableHelper&) = delete;

  protected:

    G4PhysicsTableHelper();
    ~G4PhysicsTableHelper();
 
    static G4int verboseLevel;
      // Control flag for output message
};

#endif
