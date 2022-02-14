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
// G4WorkerThread
//
// Class description:
//
// This is a class to encapsulate thread-specific data.
// Used by G4MTRunManager and G4WorkerRunManager classes.

// Authors: X.Dong, A.Dotti, 2013
// --------------------------------------------------------------------
#ifndef G4WorkerThread_hh
#define G4WorkerThread_hh 1

#include "G4String.hh"
#include "G4Types.hh"

#include "G4Threading.hh"

class G4WorkerThread
{
  public:

    void SetThreadId(G4int threadId);
    G4int GetThreadId() const;

    void SetNumberThreads(G4int numnberThreads);
    G4int GetNumberThreads() const;

    static void BuildGeometryAndPhysicsVector();
      // Build geometry for workers
    static void DestroyGeometryAndPhysicsVector();
    static void UpdateGeometryAndPhysicsVectorFromMaster();

    void SetPinAffinity(G4int aff) const;
      // Setting Pin Affinity

  private:

    G4int threadId = 0;
    G4int numThreads = 0;
};

#endif
