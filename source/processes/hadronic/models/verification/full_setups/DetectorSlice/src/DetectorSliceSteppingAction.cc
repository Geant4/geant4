#include "DetectorSliceSteppingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"


DetectorSliceSteppingAction::DetectorSliceSteppingAction() : 
  totalEdepAllParticlesInTracker( 0.0 ),
  totalEdepAllParticlesInEmCal( 0.0 ),
  totalEdepAllParticlesInHadCal( 0.0 ),
  totalEdepAllParticlesInMuonDetector( 0.0 ),
  exitingRadiusPrimaryMuon( 0.0 ),
  primaryParticleId( 0 ), primaryParticleEnergy( 0.0 ),
  isFirstStepOfTheEvent( true ),
  isPrimaryMuonReachingMuonDetector( false ),
  muonMinusId( G4MuonMinus::MuonMinusDefinition()->GetPDGEncoding() ),
  muonPlusId( G4MuonPlus::MuonPlusDefinition()->GetPDGEncoding() )
 {}


DetectorSliceSteppingAction::~DetectorSliceSteppingAction() {}


void DetectorSliceSteppingAction::UserSteppingAction( const G4Step * theStep ) {

  // Store the information about the ID and the kinetic energy 
  // of the primary particle, at the first step of the first event.
  // Notice that for the kinetic energy, we are considering the 
  // "pre-step" point of such first step.
  if ( isFirstStepOfTheEvent ) {
    if ( theStep->GetTrack()->GetParentID() == 0 ) {
      primaryParticleId = theStep->GetTrack()->GetDefinition()->GetPDGEncoding();
      primaryParticleEnergy = theStep->GetPreStepPoint()->GetKineticEnergy();
      isFirstStepOfTheEvent = false;
    }
  }

  G4double edep = theStep->GetTotalEnergyDeposit();
  if ( edep > 0.001*eV ) {
    G4String nameVolume = theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
    if ( nameVolume == "physiTracker" ) {
      totalEdepAllParticlesInTracker += edep;
    } else if ( nameVolume.find( "physiEm" ) != std::string::npos ) {
      totalEdepAllParticlesInEmCal += edep;
    } else if ( nameVolume.find( "physiHad" ) != std::string::npos ) {
      totalEdepAllParticlesInHadCal += edep;
    } else if ( nameVolume == "physiMuon" ) {
      totalEdepAllParticlesInMuonDetector += edep;
    }
  }

  // In the case of primary muon, leaving the muon detector,
  // store the exiting radius.
  if ( theStep->GetTrack()->GetParentID() == 0 ) {
    G4int particleId = theStep->GetTrack()->GetDefinition()->GetPDGEncoding();
    if ( particleId == muonMinusId  ||  particleId == muonPlusId ) {
      G4String nameVolume = theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
      if ( nameVolume == "physiMuon" ) {
	isPrimaryMuonReachingMuonDetector = true;
      } else if ( isPrimaryMuonReachingMuonDetector  &&  nameVolume == "expHall" ) {
	 exitingRadiusPrimaryMuon = 
	   std::sqrt( theStep->GetPostStepPoint()->GetPosition().x() *
		      theStep->GetPostStepPoint()->GetPosition().x() +
		      theStep->GetPostStepPoint()->GetPosition().y() *
		      theStep->GetPostStepPoint()->GetPosition().y() );  
	 isPrimaryMuonReachingMuonDetector = false;
      }
    }
  }

}


void DetectorSliceSteppingAction::reset() {
  totalEdepAllParticlesInTracker = 0.0;
  totalEdepAllParticlesInEmCal = 0.0;
  totalEdepAllParticlesInHadCal = 0.0;
  totalEdepAllParticlesInMuonDetector = 0.0;
  exitingRadiusPrimaryMuon = 0.0;
  primaryParticleId = 0;
  primaryParticleEnergy = 0.0;
  isFirstStepOfTheEvent = true;
  isPrimaryMuonReachingMuonDetector = false;
}

