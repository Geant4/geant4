// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EvtBiasMechanism.cc,v 1.3 1999-12-15 14:53:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD Group
//	
//	
// ------------------------------------------------------------
//   Implemented for the new scheme                 9 Nov. 1998  H.Kurahige
// --------------------------------------------------------------

#include "G4EvtBiasMechanism.hh"
#include "G4VParticleChange.hh"
#include "G4Track.hh"
#include "G4Step.hh"


G4EvtBiasMechanism::G4EvtBiasMechanism(const G4String& name, G4int mulFactor):
   G4VEvtBiasMechanism(name),
   MultiplicationForSecondaries(mulFactor),
   particleToBeBiased(0)
{
}  

G4EvtBiasMechanism::G4EvtBiasMechanism(const G4EvtBiasMechanism& right):
   G4VEvtBiasMechanism(right),
   MultiplicationForSecondaries(right.MultiplicationForSecondaries)
{
   particleToBeBiased = right.particleToBeBiased;
}

G4EvtBiasMechanism::~G4EvtBiasMechanism()
{
}

G4VParticleChange* G4EvtBiasMechanism::ApplyMath( G4VParticleChange* pChange, 
						  const G4Step& aStep )
{  
  if (particleToBeBiased != 0) {
    G4int currentNumberOfSecondaries = pChange->GetNumberOfSecondaries();
    G4int totalNumberOfSecondaries = currentNumberOfSecondaries;
    G4int idx;
    G4Track* track;
    G4double theParentWeight = pChange->GetParentWeight();

    G4TrackFastVector* tempList = new G4TrackFastVector();
    tempList->Initialize(currentNumberOfSecondaries);

    // fill tempList
    for (idx=0; idx<currentNumberOfSecondaries; idx+=1){
      track = pChange->GetSecondary(idx);
      tempList->SetElement(idx, track);
      if ( particleToBeBiased == track->GetDefinition() ) {
	  totalNumberOfSecondaries += (MultiplicationForSecondaries-1);
      }
    }

    pChange->Clear();
    pChange->SetNumberOfSecondaries(totalNumberOfSecondaries);
    for (idx=0; idx<currentNumberOfSecondaries; idx+=1){
      track = (*tempList)[idx];
      if (particleToBeBiased == track->GetDefinition()) {
        track->SetWeight( theParentWeight/double( MultiplicationForSecondaries) );
        pChange->AddSecondary(track);
        for (G4int i=0; i<MultiplicationForSecondaries-1; i+=1) {
          // duplicate track
          G4Track* newTrack = new G4Track(
			       new G4DynamicParticle( *(track->GetDynamicParticle()) ),
			       track->GetGlobalTime(),
			       track->GetPosition()
                             );
          newTrack->SetWeight( theParentWeight/double( MultiplicationForSecondaries) );
	  pChange->AddSecondary(newTrack);
        }
      } else {
        track->SetWeight( theParentWeight );
        pChange->AddSecondary(track);
      }
    }
    delete tempList;
  }
  return pChange;
}

