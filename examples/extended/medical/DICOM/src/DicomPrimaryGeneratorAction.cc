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
#include "DicomGeometry.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

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
    = particleTable->FindParticle(particleName="gamma");
  particleGun->SetParticleDefinition(particle);
  // ---- MGP ---- Numbers in the code should be replaced by const
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.*cm,0.*cm));       
  particleGun->SetParticleEnergy(5.*MeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  particleGun->GeneratePrimaryVertex(anEvent);
}

