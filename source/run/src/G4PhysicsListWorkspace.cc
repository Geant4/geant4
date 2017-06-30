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
#include "G4PhysicsListWorkspace.hh"

namespace
{
    G4PhysicsListWorkspace::pool_type thePool;
}

G4PhysicsListWorkspace::pool_type*
G4PhysicsListWorkspace::GetPool() { return &thePool; }

G4PhysicsListWorkspace::G4PhysicsListWorkspace(G4bool verbose)
   : fVerbose(verbose)
{
 
  fpVUPLSIM =
    &const_cast<G4VUPLManager&>(G4VUserPhysicsList::GetSubInstanceManager());
  fpVPCSIM =
    &const_cast<G4VPCManager&>(G4VPhysicsConstructor::GetSubInstanceManager());
  fpVMPLSIM =
    &const_cast<G4VMPLManager&>(G4VModularPhysicsList::GetSubInstanceManager());
    
  // Copy information from master into PolyCone/Gon Sides in this thread.
  InitialiseWorkspace();

  // Capture its address of ParticleDefinition split-class in this thread
  fpVUPLOffset = fpVUPLSIM->GetOffset();
  fpVPCOffset  = fpVPCSIM->GetOffset();
  fpVMPLOffset = fpVMPLSIM->GetOffset();
}

G4PhysicsListWorkspace::~G4PhysicsListWorkspace()
{
}

void
G4PhysicsListWorkspace::UseWorkspace()
{
    if( fVerbose )
        G4cout << "G4PhysicsListWorkspace::UseWorkspace: "
               << "Copying particles-definition Split-Class - Start " << G4endl;

    // Implementation copied from
    // G4WorkerThread::BuildGeometryAndPhysicsVector()
  
    // Physics List related, split classes mechanism:
    // instantiate sub-instance for this thread
    fpVUPLSIM->UseWorkArea(fpVUPLOffset);
    fpVPCSIM->UseWorkArea(fpVPCOffset);
    fpVMPLSIM->UseWorkArea(fpVMPLOffset);
}

void G4PhysicsListWorkspace::ReleaseWorkspace()
{
    fpVUPLSIM->UseWorkArea(0);
    fpVPCSIM->UseWorkArea(0);
    fpVMPLSIM->UseWorkArea(0);
}

void G4PhysicsListWorkspace::InitialisePhysicsList()
{
}

void
G4PhysicsListWorkspace::InitialiseWorkspace()
{
    if( fVerbose )
      G4cout << "G4PhysicsListWorkspace::InitialiseWorkspace: "
             << "Copying particles-definition Split-Class - Start " << G4endl;
    
    // PhysicsList related, split classes mechanism:
    //   Do *NOT* instantiate sub-instance for this thread,
    //   just copy the contents !!
    fpVUPLSIM->NewSubInstances();
    fpVPCSIM->NewSubInstances();
    // The following line is fundamental! If we call NewSubInstances it will not work
    // See: https://jira-geant4.kek.jp/browse/DEV-284
    fpVMPLSIM->WorkerCopySubInstanceArray();

    // Additional initialization if needed - beyond copying memory
    InitialisePhysicsList();
  
    if( fVerbose )
      G4cout << "G4PhysicsListWorkspace::CreateAndUseWorkspace: "
             << "Copying particles-definition Split-Class - Done!" << G4endl;
}

void G4PhysicsListWorkspace::DestroyWorkspace()
{
    fpVUPLSIM->FreeWorker();
    fpVPCSIM->FreeWorker();
    fpVMPLSIM->FreeWorker();
}
