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
#include "G4ParticlesWorkspace.hh"

namespace
{
  G4ParticlesWorkspace::pool_type thePool;
}

G4ParticlesWorkspace::pool_type*
G4ParticlesWorkspace::GetPool() { return &thePool; }

G4ParticlesWorkspace::G4ParticlesWorkspace(G4bool verbose)
   : fVerbose(verbose)
{
  fpPDefSIM = 
    &const_cast<G4PDefManager&>(G4ParticleDefinition::GetSubInstanceManager());

  // Copy information from master into PolyCone/Gon Sides in this thread.
  InitialiseWorkspace();

  // Capture its address of ParticleDefinition split-class in this thread
  fpPDefOffset = fpPDefSIM->GetOffset();   
}

G4ParticlesWorkspace::~G4ParticlesWorkspace()
{
}

void
G4ParticlesWorkspace::UseWorkspace()
{
  if( fVerbose ) 
     G4cout << "G4ParticlesWorkspace::UseWorkspace: "
            << "Copying particles-definition Split-Class - Start " << G4endl;

  // Implementation copied from G4WorkerThread::BuildGeometryAndPhysicsVector()
  
  // Geometry related, split classes mechanism: instantiate
  // sub-instance for this thread
  fpPDefSIM->UseWorkArea(fpPDefOffset);
}

void G4ParticlesWorkspace::ReleaseWorkspace()
{
  fpPDefSIM->UseWorkArea(0);
}

void G4ParticlesWorkspace::InitialiseParticles()
{
}

void
G4ParticlesWorkspace::InitialiseWorkspace()
{
  if( fVerbose ) 
     G4cout << "G4ParticlesWorkspace::InitialiseWorkspace: "
            << "Copying particles-definition Split-Class - Start " << G4endl;
    
  // Particles related, split classes mechanism:
  //   Do *NOT* instantiate sub-instance for this thread,
  //   just copy the contents !!
  fpPDefSIM->NewSubInstances();

  // Additional initialization if needed - beyond copying memory
  InitialiseParticles();
  
  if( fVerbose ) 
     G4cout << "G4ParticlesWorkspace::CreateAndUseWorkspace: "
            << "Copying particles-definition Split-Class - Done!" << G4endl;
}

void G4ParticlesWorkspace::DestroyWorkspace()
{
  fpPDefSIM->FreeSlave();
}
