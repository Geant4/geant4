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
//
// $Id: PhotInPrimaryGeneratorAction.cc,v 1.6 2006/06/29 16:25:21 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//

#define debug

#include "PhotInPrimaryGeneratorAction.hh"

PhotInPrimaryGeneratorAction::PhotInPrimaryGeneratorAction():
  section(1),detector(0),part("gamma"),energy(100.),oldPart("gamma"),oldEnergy(100.)
{
#ifdef debug
  G4cout<<"PhotInPrimaryGeneratorAction::Constructor: part="<<part<<", E="<<energy<<G4endl;
#endif
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle); // Initialization of the pointer from BODY

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName=part);
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(energy*GeV);
}

PhotInPrimaryGeneratorAction::~PhotInPrimaryGeneratorAction() {delete particleGun;}

void PhotInPrimaryGeneratorAction::SetProjectileName(G4String partName)
{
#ifdef debug
  G4cout<<"PhotInPrimaryGeneratorAction::SetProjectileName:Before Name="<<partName
								<<", oldName="<<part<<G4endl;
#endif
  //@@ Make a check that such a particle exists (M.K.)
  part=partName;
#ifdef debug
  G4cout<<"PhotInPrimaryGeneratorAction::SetProjectileName: After Name="<<part<<G4endl;
#endif
}

void PhotInPrimaryGeneratorAction::SetProjectileEnergy(G4double partEnergy)
{
#ifdef debug
  G4cout<<"PhotInPrimaryGeneratorAction::SetProjectileEnergy: before E="<<energy<<G4endl;
#endif
  energy=partEnergy;
#ifdef debug
  G4cout<<"PhotInPrimaryGeneratorAction::SetProjectileEnergy: before E="<<energy<<G4endl;
#endif
}
 
void PhotInPrimaryGeneratorAction::SetDetector(PhotInDetectorConstruction* det)
{
#ifdef debug
  G4cout<<"PhotInPrimaryGeneratorAction::SetDetector: is called"<<G4endl;
#endif
  detector=det;
}

void PhotInPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
#ifdef debug
  G4cout<<"PhotInPrimaryGeneratorAction::GeneratePrimaries: "<<part<<",E="<<energy<<G4endl;
#endif
  if(part!=oldPart || energy!=oldEnergy)
  {
#ifdef debug
    G4cout<<"PhotInPrimaryGeneratorAction::GeneratePrimaries: part&E're redefined"<<G4endl;
#endif
    // new particle kinematic
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle = particleTable->FindParticle(particleName=part);
    particleGun->SetParticleDefinition(particle);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    particleGun->SetParticleEnergy(energy*GeV);
    oldPart=part;
    oldEnergy=energy;
  }

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

