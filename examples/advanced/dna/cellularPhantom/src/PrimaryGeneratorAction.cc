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
// -----------------------------------------------------------------------------
//       MONTE CARLO SIMULATION OF REALISTIC GEOMETRY FROM MICROSCOPES IMAGES
//
// Authors and contributors:
// P. Barberet (a), S. Incerti (a), N. H. Tran (a), L. Morelli (a,b)
//
// a) University of Bordeaux, CNRS, LP2i, UMR5797, Gradignan, France
// b) Politecnico di Milano, Italy
//
// If you use this code, please cite the following publication:
// P. Barberet et al.,
// "Monte-Carlo dosimetry on a realistic cell monolayer
// geometry exposed to alpha particles."
// Ph. Barberet et al 2012 Phys. Med. Biol. 57 2189
// doi: 110.1088/0031-9155/57/8/2189
// -----------------------------------------------------------------------------

#include "PrimaryGeneratorAction.hh"

#include <G4GeneralParticleSource.hh>
#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
:G4VUserPrimaryGeneratorAction()
{
  fGPS = new G4GeneralParticleSource();

  G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("proton");

  fGPS->SetParticleDefinition(particle);
  fGPS->GetCurrentSource()->GetEneDist()->SetMonoEnergy(6 * MeV);
  fGPS->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  fGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0., 0., -1. * mm));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fGPS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fGPS->GeneratePrimaryVertex(anEvent);
}
