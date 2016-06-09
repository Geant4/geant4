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
// $Id: TiaraStackingAction.cc,v 1.4 2003/06/25 09:13:13 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

