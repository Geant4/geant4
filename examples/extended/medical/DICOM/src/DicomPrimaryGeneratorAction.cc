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
// The code was written by :
//	*Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//
// + Université Laval, Québec (QC) Canada
//*******************************************************

#include "DicomPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "RegularDicomDetectorConstruction.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "CLHEP/Random/RandFlat.h"

DicomPrimaryGeneratorAction::DicomPrimaryGeneratorAction()
{
  G4int nParticle = 1;
  particleGun  = new G4ParticleGun(nParticle);  	     
}

DicomPrimaryGeneratorAction::~DicomPrimaryGeneratorAction()
{
  delete particleGun;
}

void DicomPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  // put the e- in the x direction of the patient (z in the accelerator axs) to hit patient in the central slice of the phantom
  G4ThreeVector dir(0,0,1);
  //G4ThreeVector dir(2.*CLHEP::RandFlat::shoot()-1.,2.*CLHEP::RandFlat::shoot()-1.,-CLHEP::RandFlat::shoot());
  dir /= dir.mag();
  particleGun->SetParticleMomentumDirection(dir);       
  particleGun->SetParticleEnergy(10.*MeV);
  //put it at SAD = 1 m on xy plane of central slice
  particleGun->SetParticlePosition(G4ThreeVector(0.,0.,-0.1));
  //particleGun->SetParticlePosition(G4ThreeVector(0.,0.,-22.));
  particleGun->GeneratePrimaryVertex(anEvent);
}

