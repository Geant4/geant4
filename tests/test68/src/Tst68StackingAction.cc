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
#include <iostream>
#include <string>
#include <cmath>

#include "Tst68StackingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

Tst68StackingAction::Tst68StackingAction() :

  numberGammas1( 0 ), numberElectrons1( 0 ), numberPositrons1( 0 ),
  numberGammas2( 0 ), numberElectrons2( 0 ), numberPositrons2( 0 ),
  numberNeutrons( 0 ),

  gammaEnergyThreshold1( 100.0*MeV ),
  electronEnergyThreshold1( 50.0*MeV ),
  positronEnergyThreshold1( 50.0*MeV ),
  gammaWeightThreshold1( 100.0 ),
  electronWeightThreshold1( 100.0 ),
  positronWeightThreshold1( 100.0 ),
  gammaBiasingFactor1( 4 ),
  electronBiasingFactor1( 4 ),
  positronBiasingFactor1( 4 ),

  gammaEnergyThreshold2( 200.0*MeV ),
  electronEnergyThreshold2( 200.0*MeV ),
  positronEnergyThreshold2( 200.0*MeV ),
  gammaWeightThreshold2( 100.0 ),
  electronWeightThreshold2( 100.0 ),
  positronWeightThreshold2( 100.0 ),
  gammaBiasingFactor2( 2 ),
  electronBiasingFactor2( 2 ),
  positronBiasingFactor2( 2 ),

  //neutronKillingEnergyThreshold( 1.0*MeV ),
  neutronKillingEnergyThreshold( 0.0 ),
  neutronEnergyThreshold( 80.0*MeV ),
  neutronWeightThreshold( 100.0 ),
  neutronBiasingFactor( 4 )

{
  for ( int i = 0; i < Nmax; i++ ) {
    weightGammaVec[ i ] = 0;  
    weightElectronVec[ i ] = 0;  
    weightPositronVec[ i ] = 0;  
    weightNeutronVec[ i ] = 0;  
  }

  G4cout << " --- Tst68StackingAction constructor --- " << G4endl
         << " \t gammaEnergyThreshold1    = " << gammaEnergyThreshold1 / MeV 
	 << " MeV " << G4endl  
         << " \t electronEnergyThreshold1 = " << electronEnergyThreshold1 / MeV 
	 << " MeV " << G4endl  
         << " \t positronEnergyThreshold1 = " << positronEnergyThreshold1 / MeV 
	 << " MeV " << G4endl  
         << " \t gammaWeightThreshold1    = " << gammaWeightThreshold1 << G4endl  
         << " \t electronWeightThreshold1 = " << electronWeightThreshold1 << G4endl  
         << " \t positronWeightThreshold1 = " << positronWeightThreshold1 << G4endl  
         << " \t gammaBiasingFactor1      = " << gammaBiasingFactor1 << G4endl  
         << " \t electronBiasingFactor1   = " << electronBiasingFactor1 << G4endl  
         << " \t positronBiasingFactor1   = " << positronBiasingFactor1 << G4endl  
         << G4endl
         << " \t gammaEnergyThreshold2    = " << gammaEnergyThreshold2 / MeV 
	 << " MeV " << G4endl  
         << " \t electronEnergyThreshold2 = " << electronEnergyThreshold2 / MeV 
	 << " MeV " << G4endl  
         << " \t positronEnergyThreshold2 = " << positronEnergyThreshold2 / MeV 
	 << " MeV " << G4endl  
         << " \t gammaWeightThreshold2    = " << gammaWeightThreshold2 << G4endl  
         << " \t electronWeightThreshold2 = " << electronWeightThreshold2 << G4endl  
         << " \t positronWeightThreshold2 = " << positronWeightThreshold2 << G4endl  
         << " \t gammaBiasingFactor2      = " << gammaBiasingFactor2 << G4endl  
         << " \t electronBiasingFactor2   = " << electronBiasingFactor2 << G4endl  
         << " \t positronBiasingFactor2   = " << positronBiasingFactor2 << G4endl  
         << G4endl
         << " \t neutronKillingEnergyThreshold = " << neutronKillingEnergyThreshold / MeV
         << " MeV " << G4endl  
         << " \t neutronEnergyThreshold        = " << neutronEnergyThreshold / MeV
         << " MeV " << G4endl  
         << " \t neutronWeightThreshold        = " << neutronWeightThreshold << G4endl  
         << " \t neutronBiasingFactor          = " << neutronBiasingFactor << G4endl  
         << " ----------------------------------------------- " << G4endl;

  if ( gammaEnergyThreshold1 - gammaEnergyThreshold2 > 0.001*eV  ||
       gammaBiasingFactor1 < gammaBiasingFactor2 ) {
    G4cout << "\t INCONSISTENT INPUT PARAMETERS for Gammas!" << G4endl;
  }
  if ( electronEnergyThreshold1 - electronEnergyThreshold2 > 0.001*eV  ||
       electronBiasingFactor1 < electronBiasingFactor2 ) {
    G4cout << "\t INCONSISTENT INPUT PARAMETERS for Electrons!" << G4endl;
  }
  if ( positronEnergyThreshold1 - positronEnergyThreshold2 > 0.001*eV  ||
       positronBiasingFactor1 < positronBiasingFactor2 ) {
    G4cout << "\t INCONSISTENT INPUT PARAMETERS for Positrons!" << G4endl;
  }

}


Tst68StackingAction::~Tst68StackingAction() {

  std::cout << " ~Tst68StackingAction : WEIGHT DISTRIBUTION " << std::endl;

  std::cout << "\t GAMMA " << std::endl;
  G4double totNumTracks = 0;
  for ( int i = 0; i < Nmax; i++ ) {
    totNumTracks += weightGammaVec[ i ];
    std::cout << "\t exponent=" << i 
              << "\t weight=" << std::pow( static_cast<double>( gammaBiasingFactor2 ), 
					   i )
	      << "\t number of tracks=" << weightGammaVec[ i ] << std::endl;  
  }
  std::cout << " Total number of Gammas = " << totNumTracks << std::endl;

  std::cout << "\t ELECTRON " << std::endl;
  totNumTracks = 0;
  for ( int i = 0; i < Nmax; i++ ) {
    totNumTracks += weightElectronVec[ i ];
    std::cout << "\t exponent=" << i 
              << "\t weight=" << std::pow( static_cast<double>( electronBiasingFactor2 ),
					   i )
	      << "\t number of tracks=" << weightElectronVec[ i ] << std::endl;  
  }
  std::cout << " Total number of Electrons = " << totNumTracks << std::endl;

  std::cout << "\t POSITRON " << std::endl;
  totNumTracks = 0;
  for ( int i = 0; i < Nmax; i++ ) {
    totNumTracks += weightPositronVec[ i ];
    std::cout << "\t exponent=" << i 
              << "\t weight=" << std::pow( static_cast<double>( positronBiasingFactor2 ),
					   i )
	      << "\t number of tracks=" << weightPositronVec[ i ] << std::endl;  
  }
  std::cout << " Total number of Positron = " << totNumTracks << std::endl;

  std::cout << "\t NEUTRON " << std::endl;
  totNumTracks = 0;
  for ( int i = 0; i < Nmax; i++ ) {
    totNumTracks += weightNeutronVec[ i ];
    std::cout << "\t exponent=" << i 
              << "\t weight=" << std::pow( static_cast<double>( neutronBiasingFactor ),
					   i )
	      << "\t number of tracks=" << weightNeutronVec[ i ] << std::endl;  
  }
  std::cout << " Total number of Neutron = " << totNumTracks << std::endl;

}


G4ClassificationOfNewTrack 
Tst68StackingAction::ClassifyNewTrack( const G4Track * aTrack ) {

  G4ClassificationOfNewTrack result( fUrgent );

  // Because the method input parameter, aTrack, is defined as a pointer
  // to a const track, it is not possible to set the weight of this
  // track using such pointer aTrack. We need therefore to force a 
  // const_cast in order to be able to set the weight of the track
  // (this is essential for implementing the biasing).
  G4Track * theTrack = const_cast< G4Track * > (aTrack); 

  G4ParticleDefinition * particleType = aTrack->GetDefinition();

  if ( particleType == G4Gamma::GammaDefinition() ) {
    if ( aTrack->GetKineticEnergy() < gammaEnergyThreshold1  &&
         aTrack->GetWeight() < gammaWeightThreshold1 ) {
      numberGammas1++;
      if ( numberGammas1 % gammaBiasingFactor1 == 1 ) {
	theTrack->SetWeight( aTrack->GetWeight() * gammaBiasingFactor1 );
      } else {
	result = fKill;
      }
    } else if ( aTrack->GetKineticEnergy() < gammaEnergyThreshold2  &&
		aTrack->GetWeight() < gammaWeightThreshold2 ) {
      numberGammas2++;
      if ( numberGammas2 % gammaBiasingFactor2 == 1 ) {
	theTrack->SetWeight( aTrack->GetWeight() * gammaBiasingFactor2 );
      } else {
	result = fKill;
      }
    }

  } else if ( particleType == G4Electron::ElectronDefinition()  ) {
    if ( aTrack->GetKineticEnergy() < electronEnergyThreshold1  &&
	 aTrack->GetWeight() < electronWeightThreshold1 ) {
      numberElectrons1++;
      if ( numberElectrons1 % electronBiasingFactor1 == 1 ) {
	theTrack->SetWeight( aTrack->GetWeight() * electronBiasingFactor1 );
      } else {
	result = fKill;
      }
    } else if ( aTrack->GetKineticEnergy() < electronEnergyThreshold2  &&
		aTrack->GetWeight() < electronWeightThreshold2 ) {
      numberElectrons2++;
      if ( numberElectrons2 % electronBiasingFactor2 == 1 ) {
	theTrack->SetWeight( aTrack->GetWeight() * electronBiasingFactor2 );
      } else {
	result = fKill;
      }
    }

  } else if ( particleType == G4Positron::PositronDefinition() ) {
    if ( aTrack->GetKineticEnergy() < positronEnergyThreshold1  &&
	 aTrack->GetWeight() < positronWeightThreshold1 ) {
      numberPositrons1++;
      if ( numberPositrons1 % positronBiasingFactor1 == 1 ) {
	theTrack->SetWeight( aTrack->GetWeight() * positronBiasingFactor1 );
      } else {
	result = fKill;
      }
    } else if ( aTrack->GetKineticEnergy() < positronEnergyThreshold2  &&
		aTrack->GetWeight() < positronWeightThreshold2 ) {
      numberPositrons2++;
      if ( numberPositrons2 % positronBiasingFactor2 == 1 ) {
	theTrack->SetWeight( aTrack->GetWeight() * positronBiasingFactor2 );
      } else {
	result = fKill;
      }
    }

  }

  if ( particleType == G4Gamma::GammaDefinition() ) {
    G4int intWeight = static_cast< int >( aTrack->GetWeight() );
    G4int exponent = 0;
    // Rather than using the logarithm in base  gammaBiasingFactor 
    // to get the exponent, we divide iteratively by  gammaBiasingFactor
    // until it is greater than  gammaBiasingFactor : this is 
    // slower but it avoids round-off problems.
    while ( intWeight >= gammaBiasingFactor2 ) {
      exponent++;
      intWeight = intWeight / gammaBiasingFactor2;
    }
    if ( exponent > Nmax - 1 ) exponent = Nmax - 1;
    weightGammaVec[ exponent ]++;
  }

  if ( particleType == G4Electron::ElectronDefinition() ) {
    G4int intWeight = static_cast< int >( aTrack->GetWeight() );
    G4int exponent = 0;
    while ( intWeight >= electronBiasingFactor2 ) {
      exponent++;
      intWeight = intWeight / electronBiasingFactor2;
    }
    if ( exponent > Nmax - 1 ) exponent = Nmax - 1;
    weightElectronVec[ exponent ]++;
  }

  if ( particleType == G4Positron::PositronDefinition() ) {
    G4int intWeight = static_cast< int >( aTrack->GetWeight() );
    G4int exponent = 0;
    while ( intWeight >= positronBiasingFactor2 ) {
      exponent++;
      intWeight = intWeight / positronBiasingFactor2;
    }
    if ( exponent > Nmax - 1 ) exponent = Nmax - 1;
    weightPositronVec[ exponent ]++;
  }

  // Neutrons.
  if ( particleType == G4Neutron::NeutronDefinition() ) {
    if ( aTrack->GetKineticEnergy() < neutronKillingEnergyThreshold ) {
      // Simply kill them if below an energy threshold.
      result = fKill;
    } else if ( aTrack->GetKineticEnergy() < neutronEnergyThreshold  &&
                aTrack->GetWeight() < neutronWeightThreshold ) {
      numberNeutrons++;
      if ( numberNeutrons % neutronBiasingFactor == 1 ) {
	theTrack->SetWeight( aTrack->GetWeight() * neutronBiasingFactor );
      } else {
	result = fKill;
      }
    }
    G4int intWeight = static_cast< int >( aTrack->GetWeight() );
    G4int exponent = 0;
    while ( intWeight >= neutronBiasingFactor ) {
      exponent++;
      intWeight = intWeight / neutronBiasingFactor;
    }
    if ( exponent > Nmax - 1 ) exponent = Nmax - 1;
    weightNeutronVec[ exponent ]++;
  }

  return result;

}


void Tst68StackingAction::PrepareNewEvent() { 
  numberGammas1 = 0;
  numberElectrons1 = 0;
  numberPositrons1 = 0;
  numberGammas2 = 0;
  numberElectrons2 = 0;
  numberPositrons2 = 0;
  numberNeutrons = 0;
}


