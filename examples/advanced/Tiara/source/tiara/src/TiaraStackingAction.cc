// $Id: TiaraStackingAction.cc,v 1.3 2003-06-18 16:40:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

