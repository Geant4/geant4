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
// $Id: TiaraStackingAction.cc,v 1.5 2006/06/29 15:45:34 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//

#include "TiaraStackingAction.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4Track.hh"


TiaraStackingAction::
TiaraStackingAction() :
  fParticleCut(),
  fMinEnergyCut(-1)
{}

TiaraStackingAction::~TiaraStackingAction()
{}

void TiaraStackingAction::AddParticleCut(const G4String &particle, 
					 G4double cut){
  fParticleCut[particle] = cut;
  if (cut < fMinEnergyCut || fMinEnergyCut < 0) {
    fMinEnergyCut = cut;
  }
}

G4ClassificationOfNewTrack
TiaraStackingAction::ClassifyNewTrack(const G4Track* aTrack) {
  G4ClassificationOfNewTrack result(fUrgent);

  G4double kinEnergy = aTrack->GetKineticEnergy();

  if (fMinEnergyCut>0 && kinEnergy < fMinEnergyCut ) {
    // kill any particle with kinetic enerfy lower than the lowest cut
    result = fKill;
  }
  else {
    G4ParticleDefinition* pParticle = aTrack->GetDefinition();
    G4String pName = pParticle->GetParticleName();
    std::map<G4String, G4double>::const_iterator itP = 
      fParticleCut.find(pName);
    if (itP != fParticleCut.end()) {
      if (kinEnergy < itP->second) {
	// kill particle if it's energy is below it's cut
	result = fKill;
      }
    }
    else {
      // kill any particle not in the map
      result = fKill;
    }
  }
  return result;
}

