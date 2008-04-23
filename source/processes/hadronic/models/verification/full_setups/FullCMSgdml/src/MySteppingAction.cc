#include "MySteppingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"


MySteppingAction::MySteppingAction() : totalEdepAllParticles( 0.0 ) {}


MySteppingAction::~MySteppingAction() {}


void MySteppingAction::UserSteppingAction( const G4Step * theStep ) {

  // Consider the total energy depositions for any
  // particle anywhere in the whole experimental hall.
  // Notice that for simplicity we are not considering 
  // the energy deposit inside the sensitive parts of
  // the CMS detector, or even inside the whole CMS detector,
  // but overall in the experimental hall.

  G4double edep = theStep->GetTotalEnergyDeposit();
  if ( edep > 0.001*eV ) totalEdepAllParticles += edep;

}


