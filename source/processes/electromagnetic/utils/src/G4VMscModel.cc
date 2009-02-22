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
// $Id: G4VMscModel.cc,v 1.9 2009-02-22 17:32:08 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4VMscModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 07.03.2008
//
// Modifications:
//
//
// Class Description:
//
// General interface to msc models

// -------------------------------------------------------------------
//

#include "G4VMscModel.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4TransportationManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VMscModel::G4VMscModel(const G4String& nam):
  G4VEmModel(nam), 
  safetyHelper(0),
  facrange(0.02),
  facgeom(2.5),
  facsafety(0.25),
  skin(3.0),
  dtrl(0.05),
  lambdalimit(mm),
  geommax(1.e50*mm),
  steppingAlgorithm(fUseSafety),
  samplez(false),
  latDisplasment(true)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VMscModel::~G4VMscModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VMscModel::InitialiseSafetyHelper()
{
  if(!safetyHelper) {
    safetyHelper = G4TransportationManager::GetTransportationManager()
      ->GetSafetyHelper();
    safetyHelper->InitialiseHelper();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VMscModel::ComputeDisplacement(G4ParticleChangeForMSC* fParticleChange,  
				      const G4ThreeVector& dir,
				      G4double displacement,
				      G4double postsafety)
{
  const G4ThreeVector* pos = fParticleChange->GetProposedPosition();
  G4double r = displacement;
  if(r >  postsafety) {
    G4double newsafety = safetyHelper->ComputeSafety(*pos);
    if(r > newsafety) r = newsafety;
  }
  if(r > 0.) {

    // compute new endpoint of the Step
    G4ThreeVector newPosition = *pos + r*dir;

    // definitely not on boundary
    if(displacement == r) {
      safetyHelper->ReLocateWithinVolume(newPosition);

    } else {
      // check safety after displacement
      G4double postsafety = safetyHelper->ComputeSafety(newPosition);

      // displacement to boundary
      if(postsafety <= 0.0) {
	safetyHelper->Locate(newPosition, 
			     *fParticleChange->GetProposedMomentumDirection());

	// not on the boundary
      } else {
	safetyHelper->ReLocateWithinVolume(newPosition);
      }
    }
    fParticleChange->ProposePosition(newPosition);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
