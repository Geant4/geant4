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
//    **************************************
//    *                                    *
//    *    CellPrimaryGeneratorAction.cc   *
//    *                                    *
//    **************************************
//
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//	   Barbara Mascialino (Barbara.Mascialino@ge.infn.it)
//
// History:
// -----------
// 20 September 2006   S. Guatelli, B. Mascialino   1st implementation
//
// -------------------------------------------------------------------

#include "G4IonTable.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "CellPrimaryGeneratorAction.hh"
#include "CellAnalysisManager.hh"

CellPrimaryGeneratorAction::CellPrimaryGeneratorAction()
{
  G4int particleNumber = 1;

  particleGun = new G4ParticleGun(particleNumber);
  
  // Definition of the default particles

  // Type of primary particle
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable -> FindParticle("e-");
  particleGun -> SetParticleDefinition(particle);

  // Direction of the primary vertex : along z axis
  particleGun -> SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
 
  // Energy of the primary particle
  particleGun -> SetParticleEnergy(1.*MeV);

  // Initial position of the primary particle
// BABS Sept 20th 2006
//  particleGun -> SetParticlePosition(G4ThreeVector(0.*m,0.*m,-1.*m));
  particleGun -> SetParticlePosition(G4ThreeVector(0.*m,0.*m,-1.*m));
}

CellPrimaryGeneratorAction::~CellPrimaryGeneratorAction()
{
  delete particleGun;
}

void CellPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // The energy of the primary particles is stored in the histogram  
   CellAnalysisManager* analysis = CellAnalysisManager::getInstance();
   analysis -> primaryparticle_energy(GetInitialEnergy());

  // Primary particle is generated
  particleGun -> GeneratePrimaryVertex(anEvent);
}

G4double CellPrimaryGeneratorAction::GetInitialEnergy()
{
  G4double primaryParticleEnergy = particleGun -> GetParticleEnergy(); 
  return primaryParticleEnergy;
}

G4String CellPrimaryGeneratorAction::GetParticle()
{
  G4String primaryParticleName = particleGun->GetParticleDefinition()->GetParticleName();
  return primaryParticleName;
}



