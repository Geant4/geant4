#include "StatAccepTestStackingAction.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include <iostream>
#include <string>


StatAccepTestStackingAction::StatAccepTestStackingAction() :
  numberGammas( 0 ), numberElectrons( 0 ), numberPositrons( 0 ),
  //***LOOKHERE***
  gammaEnergyThreshold( 100.0*MeV ),
  electronEnergyThreshold( 100.0*MeV ),
  positronEnergyThreshold( 100.0*MeV ),
  gammaWeightThreshold( 100.0 ),
  electronWeightThreshold( 100.0 ),
  positronWeightThreshold( 100.0 ),
  neutronEnergyThreshold( 1.0*MeV )
  //***endLOOKHERE***
{

  //ALB const std::string nameFileOut = "weights.txt";
  //ALB filePtr = new std::ofstream( nameFileOut.c_str(), std::ios::out ); 
  //ALB if ( ! filePtr ) {
  //ALB   std::cout << " ***ERROR*** : StatAccepTestStackingAction : CANNOT OPEN OUTPUT FILE! "
  //ALB 	      << std::endl;
  //ALB }

  for ( int i = 0; i < Nmax; i++ ) {
    weightVec[ i ] = 0;  
  }

  G4cout << " --- StatAccepTestStackingAction constructor --- " << G4endl
         << " \t gammaEnergyThreshold = "    << gammaEnergyThreshold / MeV 
	 << " MeV " << G4endl  
         << " \t electronEnergyThreshold = " << electronEnergyThreshold / MeV 
	 << " MeV " << G4endl  
         << " \t positronEnergyThreshold = " << positronEnergyThreshold / MeV 
	 << " MeV " << G4endl  
         << " \t gammaWeightThreshold = "    << gammaWeightThreshold << G4endl  
         << " \t electronWeightThreshold = " << electronWeightThreshold << G4endl  
         << " \t positronWeightThreshold = " << positronWeightThreshold << G4endl  
         << " \t neutronEnergyThreshold = " << neutronEnergyThreshold / MeV
         << " MeV " << G4endl  
         << " ----------------------------------------------- " << G4endl;

}


StatAccepTestStackingAction::~StatAccepTestStackingAction() {

  //ALB if ( filePtr ) {
  //ALB   filePtr->close();
  //ALB   delete filePtr;
  //ALB }
  
  std::cout << " ~StatAccepTestStackingAction : WEIGHT DISTRIBUTION " << std::endl;
  G4double totNumTracks = 0;
  for ( int i = 0; i < Nmax; i++ ) {
    totNumTracks += weightVec[ i ];
    std::cout << " exponent=" << i 
              << "\t weight=" << pow(2.0, i)
	      << "\t number of tracks=" << weightVec[ i ] << std::endl;  
  }
  std::cout << " Total number of tracks = " << totNumTracks << std::endl;

}


G4ClassificationOfNewTrack 
StatAccepTestStackingAction::ClassifyNewTrack(const G4Track * aTrack) {

  G4ClassificationOfNewTrack result( fUrgent );

  // Because the method input parameter, aTrack, is defined as a pointer
  // to a const track, it is not possible to set the weight of this
  // track using such pointer aTrack. We need therefore to force a 
  // const_cast in order to be able to set the weight of the track
  // (this is essential for implementing the biasing).
  G4Track * theTrack = const_cast< G4Track * > (aTrack);

  G4ParticleDefinition * particleType = aTrack->GetDefinition();

  if ( particleType == G4Gamma::GammaDefinition()         &&
       aTrack->GetKineticEnergy() < gammaEnergyThreshold  &&
       aTrack->GetWeight() < gammaWeightThreshold ) {
    numberGammas++;
    if ( numberGammas % 2 == 0 ) {
      result = fKill;
    } else {
      theTrack->SetWeight( aTrack->GetWeight() * 2.0 );
    }
  } else if ( particleType == G4Electron::ElectronDefinition()      &&
              aTrack->GetKineticEnergy() < electronEnergyThreshold  &&
              aTrack->GetWeight() < electronWeightThreshold ) {
    numberElectrons++;
    if ( numberElectrons % 2 == 0 ) {
      result = fKill;
    } else {
      theTrack->SetWeight( aTrack->GetWeight() * 2.0 );
    }
  } else if ( particleType == G4Positron::PositronDefinition()      &&
              aTrack->GetKineticEnergy() < positronEnergyThreshold  &&
              aTrack->GetWeight() < positronWeightThreshold ) {
    numberPositrons++;
    if ( numberPositrons % 2 == 0 ) {
      result = fKill;
    } else {
      theTrack->SetWeight( aTrack->GetWeight() * 2.0 );
    }
  }

  // Keep information about the weight of all tracks.
  G4int intWeight = static_cast< int >( aTrack->GetWeight() );
  G4int exponent = 0;
  // Rather than using the logarithm in base 2 to get the exponent,
  // divide iteratively by 2 until it is greater than 2 : this is 
  // slower but it avoids round-off problems.
  while ( intWeight > 1 ) {
    exponent++;
    intWeight = intWeight / 2;
  }
  if ( exponent > Nmax - 1 ) exponent = Nmax - 1;
  weightVec[ exponent ]++;
  //ALB if ( filePtr ) {
  //ALB if ( filePtr && aTrack->GetWeight() > 100.0 ) {
  //ALB    (*filePtr) << aTrack->GetWeight() << std::endl;
  //ALB }

  // Special for neutron: simply kill them if below an energy threshold.
  if ( particleType == G4Neutron::NeutronDefinition()  &&
       aTrack->GetKineticEnergy() < neutronEnergyThreshold ) {
      result = fKill;
  }

  return result;

}


void StatAccepTestStackingAction::PrepareNewEvent() { 
  numberGammas = 0;
  numberElectrons = 0;
  numberPositrons = 0;
}


