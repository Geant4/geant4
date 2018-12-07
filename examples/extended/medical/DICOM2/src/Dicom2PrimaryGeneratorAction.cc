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
/// \file Dicom2PrimaryGeneratorAction.cc
/// \brief Implementation of the Dicom2PrimaryGeneratorAction class

#include "Dicom2PrimaryGeneratorAction.hh"
#include "DicomDetectorConstruction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Dicom2PrimaryGeneratorAction::Dicom2PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(new G4ParticleGun(1)),
  fEnvelopeBox(nullptr),
  fEnvelopeVol(nullptr),
  fPosCenter(G4ThreeVector(0.0)),
  fPosDelta(G4ThreeVector(0.0)),
  fGeomFactor(0.8)
{
    // default particle kinematic
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle
            = particleTable->FindParticle(particleName="e-");
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fParticleGun->SetParticleEnergy(10.*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Dicom2PrimaryGeneratorAction::~Dicom2PrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dicom2PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    //this function is called at the begining of ecah event
    //

    // In order to avoid dependence of PrimaryGeneratorAction
    // on DetectorConstruction class we get phantomContainer volume
    // from G4LogicalVolumeStore.
    if (!fEnvelopeBox)
    {
        G4LogicalVolume* envLV
                = G4LogicalVolumeStore::GetInstance()->GetVolume("phantomContainer");
        if(envLV)
            fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
    }

    if(!fEnvelopeVol)
    {
        fEnvelopeVol = G4PhysicalVolumeStore::GetInstance()->GetVolume("phantomContainer");
    }

    // get the center of the phantom
    if(fEnvelopeVol)
    {
        fPosCenter = fEnvelopeVol->GetObjectTranslation();
    }
    else
    {
        G4ExceptionDescription msg;
        msg << "Envelope physical volume not found.\n";
        msg << "The gun will be place at the center.";
        G4Exception("Dicom2PrimaryGeneratorAction::GeneratePrimaries()",
                    "DICOM20002",JustWarning, msg);
    }

    // get the size of the phantom
    if(fEnvelopeBox)
    {
        fPosDelta = G4ThreeVector(2.0 * fEnvelopeBox->GetXHalfLength(),
                                  2.0 * fEnvelopeBox->GetYHalfLength(),
                                  2.0 * fEnvelopeBox->GetZHalfLength());
    }
    else
    {
        G4ExceptionDescription msg;
        msg << "Envelope volume of box shape not found.\n";
        msg << "The gun will be place at the center.";
        G4Exception("Dicom2PrimaryGeneratorAction::GeneratePrimaries()",
                    "DICOM20003",JustWarning, msg);
    }

    // start at center
    G4ThreeVector pos = fPosCenter;
    // get a random displacement
    G4ThreeVector delta(fPosDelta.x() * (G4UniformRand()-0.5),
                        fPosDelta.y() * (G4UniformRand()-0.5),
                        fPosDelta.z() * (G4UniformRand()-0.5));
    // scale by geom factor
    delta *= fGeomFactor;

    //G4cout << "Starting position: " << G4BestUnit(pos + delta, "Length")
    //       << ", pos = " << G4BestUnit(fPosCenter, "Length")
    //       << ", size = " << G4BestUnit(fPosDelta, "Length")
    //       << ", delta = " << G4BestUnit(delta, "Length")
    //       << G4endl;

    fParticleGun->SetParticlePosition(pos + delta);
    fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

