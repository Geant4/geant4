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
// -------------------------------------------------------------------
// $Id: MicrodosimetryPrimaryGeneratorAction.cc,v 1.1 2007-10-09 08:00:29 sincerti Exp $
// -------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"

#include "MicrodosimetryPrimaryGeneratorAction.hh"
#include "MicrodosimetryDetectorConstruction.hh"

//DNA
#include "G4DNAGenericIonsManager.hh"
//END DNA

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrodosimetryPrimaryGeneratorAction::MicrodosimetryPrimaryGeneratorAction(MicrodosimetryDetectorConstruction* DC)
  :Detector(DC)
{
  particleGun  = new G4ParticleGun(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrodosimetryPrimaryGeneratorAction::~MicrodosimetryPrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrodosimetryPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4int numEvent;
  numEvent=anEvent->GetEventID()+1;
  G4double x0,y0,z0,xMom0,yMom0,zMom0,e0;

  x0=0;
  y0=0;
  z0=0;
  xMom0=0;
  yMom0=0;
  zMom0=1;
  e0=2.37*MeV;
  e0=1.5*keV;

  // VERBOSE
  
  G4cout 
  << "-> Event # " << numEvent 
  << " generated " 
  << G4endl;
  
  particleGun->SetParticleEnergy(e0);

  particleGun->SetParticleMomentumDirection(G4ThreeVector(xMom0,yMom0,zMom0));

  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  
  G4ParticleDefinition* particle=
    G4ParticleTable::GetParticleTable()->FindParticle("alpha+");

  particleGun->SetParticleDefinition(particle);
  
  particleGun->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


