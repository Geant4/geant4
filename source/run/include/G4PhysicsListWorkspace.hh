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
// G4PhysicsListWorkspace
//
// Class description:
//
// Manage the per-thread state of lists - those which have a per-thread
// state and dependent classes (if any). In particular it
//  - owns the arrays that implement 'split' classes
//  - classes/objects which are owned by the split classes.
// The classes/objects affected are:
//  - 'split' classes part of its state is per-thread,
//  - per-thread objects, in particular those owned by the split classes.
// Goal: take ownership and control of per-thread state of classes
// to work with multi-threading.

// Authors: J.Apostolakis, A.Dotti - 4 October 2013
// --------------------------------------------------------------------
#ifndef G4PhysicsListWorkspace_hh
#define G4PhysicsListWorkspace_hh 1

#include "G4TWorkspacePool.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4VUserPhysicsList.hh"

class G4PhysicsListWorkspace
{
  public:
    using pool_type = G4TWorkspacePool<G4PhysicsListWorkspace>;

  public:
    G4PhysicsListWorkspace(G4bool verbose = false);
    ~G4PhysicsListWorkspace() = default;

    // Take ownership
    void UseWorkspace();

    // Release ownership
    void ReleaseWorkspace();

    // Release ownership and destroy
    void DestroyWorkspace();

    // To be called at start of each run (especially 2nd and further runs)
    void InitialiseWorkspace();

    inline void SetVerbose(G4bool v) { fVerbose = v; }
    inline G4bool GetVerbose() { return fVerbose; }

    static pool_type* GetPool();

  protected:  // Implementation methods
    void InitialisePhysicsList();

  private:  // Helper pointers - can be per instance or shared
    // Store SubInstanceManager object pointers (SIM pointers)
    G4VUPLManager* fpVUPLSIM = nullptr;
    G4VPCManager* fpVPCSIM = nullptr;
    G4VMPLManager* fpVMPLSIM = nullptr;

    // Per Instance variables
    // The ownership of the Data Arrays is IN this object
    G4VUPLData* fpVUPLOffset = nullptr;
    G4VPCData* fpVPCOffset = nullptr;
    G4VMPLData* fpVMPLOffset = nullptr;

    G4bool fVerbose = false;
};

#endif
