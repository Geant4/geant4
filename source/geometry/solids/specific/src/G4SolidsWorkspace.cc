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

#include "G4SolidsWorkspace.hh"

#include "G4VSolid.hh"

#include "G4PolyconeSide.hh"
#include "G4PolyhedraSide.hh"

#include "G4AutoLock.hh"

G4SolidsWorkspace::G4SolidsWorkspace(G4bool verbose)
   : fVerbose(verbose)
{
  // G4PlSideManager *fPolyconeSideSIM; // (G4PlSideManager::GetSubInstanceManager())
  
  fpPolyconeSideSIM=
      &const_cast<G4PlSideManager&>(G4PolyconeSide::GetSubInstanceManager());
  fpPolyhedraSideSIM=
      &const_cast<G4PhSideManager&>(G4PolyhedraSide::GetSubInstanceManager());

  // Copy information from master into PolyCone/Gon Sides in this thread.
  InitialiseWorkspace();

  // Capture its address of PolyCone/Gon Sides in this thread
  fPolyconeSideOffset= fpPolyconeSideSIM->GetOffset();                                
  fPolyhedraSideOffset= fpPolyhedraSideSIM->GetOffset();
}

G4SolidsWorkspace::~G4SolidsWorkspace()
{
}

//  Static methods 
//      For with current (original) G4WorkerThread -- which uses static methods

void
G4SolidsWorkspace::UseWorkspace()
{
  if( fVerbose ) 
     G4cout << "G4SolidsWorkspace::UseWorkspace: Copying geometry - Start " << G4endl;

  // Implementation copied from  G4WorkerThread::BuildGeometryAndPhysicsVector()
  
  //Geometry related, split classes mechanism: instantiate sub-instance for this thread
  fpPolyconeSideSIM->UseWorkArea(fPolyconeSideOffset);
  fpPolyhedraSideSIM->UseWorkArea(fPolyhedraSideOffset);
}


void G4SolidsWorkspace::ReleaseWorkspace()
//  The opposite of Use Workspace - let go of it.
{
  fpPolyconeSideSIM->UseWorkArea(0);
  fpPolyhedraSideSIM->UseWorkArea(0);
}

void G4SolidsWorkspace::InitialiseSolids()
{
}

void
G4SolidsWorkspace::InitialiseWorkspace()
{
  if( fVerbose ) 
     G4cout << "G4SolidsWorkspace::InitialiseWorkspace: "
            << "Copying geometry - Start " << G4endl;
    
  //Geometry related, split classes mechanism:
  //   Do *NOT* instantiate sub-instance for this thread,
  //     just copy the contents !!

  fpPolyconeSideSIM->SlaveInitializeSubInstance();
  fpPolyhedraSideSIM->SlaveInitializeSubInstance();

  // Additional initialization if needed - beyond copying memory
  InitialiseSolids();
  
  if( fVerbose ) 
     G4cout << "G4SolidsWorkspace::CreateAndUseWorkspace: "
            << "Copying geometry - Done!" << G4endl;
}

void G4SolidsWorkspace::DestroyWorkspace()
{
  fpPolyconeSideSIM->FreeSlave();
  fpPolyhedraSideSIM->FreeSlave();
}
