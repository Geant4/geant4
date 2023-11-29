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
// G4SolidsWorkspace
//
// Class description:
//
//   Manages the per-thread state of solids - those which 
//   have a per-thread state and dependent classes (if any)
//   In particular it 
//       - owns the arrays that implement 'split' classes 
//       - classes/objects which are owned by the split classes.
//   Background: the classes/objects affected are  
//       - 'split' classes part of its state is per-thread,
//       - per-thread objects, in particular those which are owned 
//         by the split classes.
// Goal: Take ownership and control of per-thread state of 
//       classes to work with multi-threading. 
//       Offshoot of G4GeometryWorkspace, to deal with Solids.

// 4.10.2013 - Created: John Apostolakis, Andrea Dotti
// --------------------------------------------------------------------
#ifndef G4SOLIDSWORKSPACE_HH
#define G4SOLIDSWORKSPACE_HH

#include "G4TWorkspacePool.hh"
#include "G4VSolid.hh"

#include "G4PolyconeSide.hh"
#include "G4PolyhedraSide.hh"

class G4SolidsWorkspace
{
  public: 

    using pool_type = G4TWorkspacePool<G4SolidsWorkspace>;

    G4SolidsWorkspace(G4bool verbose = false);
    ~G4SolidsWorkspace() = default;

    void UseWorkspace();     // Take ownership
    void ReleaseWorkspace(); // Release ownership
    void DestroyWorkspace(); // Release ownership and destroy

    void InitialiseWorkspace();
      // To be called at start of each run (especially 2nd and further runs)

    inline void SetVerbose(G4bool v) { fVerbose=v; } 
    inline G4bool GetVerbose()  { return fVerbose; } 

    static pool_type* GetPool();

  protected:

    void InitialiseSolids();

  private:

    // Helper pointers - can be per instance or shared
    //
    G4PlSideManager* fpPolyconeSideSIM = nullptr;
    G4PhSideManager* fpPolyhedraSideSIM = nullptr;
  
    // Per Instance variables
    // NOTE: the ownership of the Data Arrays is IN this object
 
    // Store SubInstanceManager object pointers (SIM pointers)
    //
    G4PlSideData* fPolyconeSideOffset = nullptr;
    G4PhSideData* fPolyhedraSideOffset = nullptr;

    G4bool fVerbose = false;
};

#endif // G4SOLIDSWORKSPACE_HH
