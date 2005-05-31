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
// $Id: PhotInPrimaryGeneratorAction.cc,v 1.2 2005-05-31 15:23:01 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#define debug

#include "PhotInPrimaryGeneratorAction.hh"

PhotInPrimaryGeneratorAction::PhotInPrimaryGeneratorAction():
  section(1),detector(0)
{
#ifdef debug
  G4cout<<"PhotInPrimaryGeneratorAction::Constructor: is called"<<G4endl;
#endif
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle); // Initialization of the pointer from BODY
  
  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="gamma");//@@
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(100.*GeV); // @@ Make a possibility to change on flight
}

PhotInPrimaryGeneratorAction::~PhotInPrimaryGeneratorAction() {delete particleGun;}

void PhotInPrimaryGeneratorAction::SetDetector(PhotInDetectorConstruction* det)
{
#ifdef debug
  G4cout<<"PhotInPrimaryGeneratorAction::SetDetector: is called"<<G4endl;
#endif
  detector=det;
}

void PhotInPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  
  G4double hY=detector->GetHalfYWidth();
  G4double hZ=detector->GetHalfZThickness();
  if(section==1) // The section #1 is always in the center (comment or change if #ofSec!=3)
  {
    particleGun->SetParticlePosition(G4ThreeVector(0.,0.,-hZ));
    particleGun->GeneratePrimaryVertex(anEvent);
  }
  else if(detector->IsSerial())
  {
    particleGun->SetParticlePosition(G4ThreeVector(0.,0.,hZ*(section+section-3)));
    particleGun->GeneratePrimaryVertex(anEvent);
  }
  else
  {
    particleGun->SetParticlePosition(G4ThreeVector(0.,hY*(section+section-2),-hZ));
    particleGun->GeneratePrimaryVertex(anEvent);
  }
}

