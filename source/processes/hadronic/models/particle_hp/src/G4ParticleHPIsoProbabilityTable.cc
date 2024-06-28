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
// -------------------------------------------------------------------
//
//      Geant4 source file 
//
//      File name: G4ParticleHPIsoProbabilityTable.cc
//
//      Authors: Marek Zmeskal (CTU, Czech Technical University in Prague, Czech Republic)
//               Loic Thulliez (CEA France)      
//
//      Creation date: 4 June 2024
//
//      Description: Class for the probability table of the given isotope
//                   and for the given temperature.
//                   It reads the files with probability tables and
//                   finds the correct cross-section.
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//

#include "G4ParticleHPIsoProbabilityTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Exception.hh"
#include "G4NucleiProperties.hh"
#include "G4ReactionProduct.hh"
#include "G4ParticleHPManager.hh"
#include "G4ParticleHPVector.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4Nucleus.hh"
#include "G4Neutron.hh"
#include "G4ParticleHPChannel.hh"
#include "G4ParticleHPChannelList.hh"
#include <string>
#include <fstream>

///--------------------------------------------------------------------------------------
G4ParticleHPIsoProbabilityTable::G4ParticleHPIsoProbabilityTable() : 
  theEnergies( nullptr ), theProbabilities( nullptr ), theElasticData( nullptr ), 
  theCaptureData( nullptr ), theFissionData( nullptr ) 
{
  Z = 0;
  A = 0;
  m = -1;
  T = -1.0;
  Emin = DBL_MAX;
  Emax = 0.0;
}

///--------------------------------------------------------------------------------------
G4ParticleHPIsoProbabilityTable::~G4ParticleHPIsoProbabilityTable() {
  for ( std::vector< std::vector< G4double >* >::iterator it = theProbabilities->begin(); 
        it != theProbabilities->end(); ++it ) {
    delete* it;
  }
  for ( std::vector< std::vector< G4double >* >::iterator it = theElasticData->begin(); 
        it != theElasticData->end(); ++it ) {
    delete* it;
  }
  for ( std::vector< std::vector< G4double >* >::iterator it = theCaptureData->begin(); 
        it != theCaptureData->end(); ++it ) {
    delete* it;
  }
  for ( std::vector< std::vector< G4double >* >::iterator it = theFissionData->begin(); 
        it != theFissionData->end(); ++it ) {
    delete* it;
  }
  delete theEnergies;
  delete theProbabilities;
  delete theElasticData;
  delete theCaptureData;
  delete theFissionData;
}

///--------------------------------------------------------------------------------------
void G4ParticleHPIsoProbabilityTable::Init( G4int, G4int, G4int, G4double, G4String ) {}

///--------------------------------------------------------------------------------------
G4double G4ParticleHPIsoProbabilityTable::GetCorrelatedIsoCrossSectionPT( const G4DynamicParticle*, G4int, 
  const G4Element*, G4double&, G4double&, std::thread::id& ) 
{
  G4Exception( "G4ParticleHPIsoProbabilityTable::GetCorrelatedIsoCrossSectionPT", "hadhp01",
	       FatalException, "The base method G4ParticleHPIsoProbabilityTable::GetCorrelatedIsoCrossSectionPT"
	       "is trying to be accessed, which is not allowed,"
	       "the derived class for NJOY or CALENDF was not properly declared." );
  return 0.0;
}

///--------------------------------------------------------------------------------------
G4double G4ParticleHPIsoProbabilityTable::GetIsoCrossSectionPT( const G4DynamicParticle*, G4int, 
  const G4Element*, G4double&, std::map< std::thread::id, G4double >&, std::thread::id& ) 
{
  G4Exception( "G4ParticleHPIsoProbabilityTable::GetIsoCrossSectionPT", "hadhp01",
	       FatalException, "The base method G4ParticleHPIsoProbabilityTable::GetIsoCrossSectionPT"
	       "is trying to be accessed, which is not allowed,"
	       "the derived class for NJOY or CALENDF was not properly declared." );
  return 0.0;
}

///--------------------------------------------------------------------------------------
G4double G4ParticleHPIsoProbabilityTable::GetDopplerBroadenedElasticXS( const G4DynamicParticle* dP, G4int indexEl, G4int isotopeJ ) {
  G4double result = 0.0;
  G4ReactionProduct theNeutron( dP->GetDefinition() );
  theNeutron.SetMomentum( dP->GetMomentum() );
  theNeutron.SetKineticEnergy( dP->GetKineticEnergy() );
  // prepare thermal nucleus
  G4Nucleus aNuc;
  G4double eleMass = G4NucleiProperties::GetNuclearMass(A, Z) / G4Neutron::Neutron()->GetPDGMass();
  G4ReactionProduct boosted;
  G4double aXsection;
  G4int counter = 0;
  G4double buffer = 0.0;
  G4int size = G4int( std::max( 10.0, T / 60.0 * kelvin ) );
  G4ThreeVector neutronVelocity = 1.0 / G4Neutron::Neutron()->GetPDGMass() * theNeutron.GetMomentum();
  G4double neutronVMag = neutronVelocity.mag();
  while ( counter == 0  ||  std::abs( buffer - result / std::max(1, counter) ) > 0.03 * buffer ) {
    if ( counter ) buffer = result / counter;
    while ( counter < size ) {
      ++counter;
      G4ReactionProduct aThermalNuc = aNuc.GetThermalNucleus( eleMass, T );
      boosted.Lorentz( theNeutron, aThermalNuc );
      G4double theEkin = boosted.GetKineticEnergy();
      aXsection = (*G4ParticleHPManager::GetInstance()->GetElasticFinalStates())[indexEl]->GetWeightedXsec( theEkin, isotopeJ );
      // velocity correction
      G4ThreeVector targetVelocity = 1.0 / aThermalNuc.GetMass() * aThermalNuc.GetMomentum();
      aXsection *= (targetVelocity - neutronVelocity).mag() / neutronVMag;
      result += aXsection;
    }
    size += size;
  }
  result /= counter;
  return result;
}

///--------------------------------------------------------------------------------------
G4double G4ParticleHPIsoProbabilityTable::GetDopplerBroadenedCaptureXS( const G4DynamicParticle* dP, G4int indexEl, G4int isotopeJ ) {
  G4double result = 0.0;
  G4ReactionProduct theNeutron( dP->GetDefinition() );
  theNeutron.SetMomentum( dP->GetMomentum() );
  theNeutron.SetKineticEnergy( dP->GetKineticEnergy() );
  // prepare thermal nucleus
  G4Nucleus aNuc;
  G4double eleMass = G4NucleiProperties::GetNuclearMass(A, Z) / G4Neutron::Neutron()->GetPDGMass();
  G4ReactionProduct boosted;
  G4double aXsection;
  G4int counter = 0;
  G4double buffer = 0.0;
  G4int size = G4int( std::max( 10.0, T / 60.0 * kelvin));
  G4ThreeVector neutronVelocity = 1.0 / G4Neutron::Neutron()->GetPDGMass() * theNeutron.GetMomentum();
  G4double neutronVMag = neutronVelocity.mag();
  while ( counter == 0  ||  std::abs( buffer - result / std::max( 1, counter ) ) > 0.03 * buffer ) {
    if ( counter ) buffer = result / counter;
    while ( counter < size ) {
      ++counter;
      G4ReactionProduct aThermalNuc = aNuc.GetThermalNucleus( eleMass, T );
      boosted.Lorentz( theNeutron, aThermalNuc );
      G4double theEkin = boosted.GetKineticEnergy();
      aXsection = (*G4ParticleHPManager::GetInstance()->GetCaptureFinalStates())[indexEl]->GetWeightedXsec( theEkin, isotopeJ );
      // velocity correction
      G4ThreeVector targetVelocity = 1.0 / aThermalNuc.GetMass() * aThermalNuc.GetMomentum();
      aXsection *= (targetVelocity - neutronVelocity).mag() / neutronVMag;
      result += aXsection;
    }
    size += size;
  }
  result /= counter;
  return result;
}

///--------------------------------------------------------------------------------------
G4double G4ParticleHPIsoProbabilityTable::GetDopplerBroadenedFissionXS( const G4DynamicParticle* dP, G4int indexEl, G4int isotopeJ ) {
  G4double result = 0.0;
  G4ReactionProduct theNeutron( dP->GetDefinition() );
  theNeutron.SetMomentum( dP->GetMomentum() );
  theNeutron.SetKineticEnergy( dP->GetKineticEnergy() );
  // prepare thermal nucleus
  G4Nucleus aNuc;
  G4double eleMass;
  eleMass = G4NucleiProperties::GetNuclearMass(A, Z) / G4Neutron::Neutron()->GetPDGMass();
  G4ReactionProduct boosted;
  G4double aXsection;
  G4int counter = 0;
  G4double buffer = 0.0;
  G4int size = G4int( std::max( 10.0, T / 60.0 * kelvin ) );
  G4ThreeVector neutronVelocity = 1.0 / G4Neutron::Neutron()->GetPDGMass() * theNeutron.GetMomentum();
  G4double neutronVMag = neutronVelocity.mag();
  while ( counter == 0  ||  std::abs( buffer - result / std::max( 1, counter ) ) > 0.01 * buffer ) {
    if ( counter ) buffer = result / counter;
    while ( counter < size ) {
      ++counter;
      G4ReactionProduct aThermalNuc = aNuc.GetThermalNucleus( eleMass, T );
      boosted.Lorentz( theNeutron, aThermalNuc );
      G4double theEkin = boosted.GetKineticEnergy();
      aXsection = (*G4ParticleHPManager::GetInstance()->GetFissionFinalStates())[indexEl]->GetWeightedXsec( theEkin, isotopeJ );
      // velocity correction
      G4ThreeVector targetVelocity = 1.0 / aThermalNuc.GetMass() * aThermalNuc.GetMomentum();
      aXsection *= (targetVelocity - neutronVelocity).mag() / neutronVMag;
      result += aXsection;
    }
    size += size;
  }
  result /= counter;
  return result;
}

///--------------------------------------------------------------------------------------
G4double G4ParticleHPIsoProbabilityTable::GetDopplerBroadenedInelasticXS( const G4DynamicParticle* dP, G4int indexEl, G4int isotopeJ ) {
  G4double result = 0.0;
  G4ReactionProduct theNeutron( dP->GetDefinition() );
  theNeutron.SetMomentum( dP->GetMomentum() );
  theNeutron.SetKineticEnergy( dP->GetKineticEnergy() );
  // prepare thermal nucleus
  G4Nucleus aNuc;
  G4double eleMass = G4NucleiProperties::GetNuclearMass(A, Z) / G4Neutron::Neutron()->GetPDGMass();
  G4ReactionProduct boosted;
  G4double aXsection;
  G4int counter = 0;
  G4int failCount = 0;
  G4double buffer = 0.0;
  G4int size = G4int( std::max( 10.0, T / 60.0 * kelvin ) );
  G4ThreeVector neutronVelocity = 1.0 / G4Neutron::Neutron()->GetPDGMass() * theNeutron.GetMomentum();
  G4double neutronVMag = neutronVelocity.mag();
  #ifndef G4PHP_DOPPLER_LOOP_ONCE
  while ( counter == 0  ||  std::abs( buffer - result / std::max( 1, counter ) ) > 0.01 * buffer ) {
    if ( counter ) buffer = result / counter;
    while ( counter < size ) {
      ++counter;
  #endif
      G4ReactionProduct aThermalNuc = aNuc.GetThermalNucleus( eleMass / G4Neutron::Neutron()->GetPDGMass(), T );
      boosted.Lorentz( theNeutron, aThermalNuc );
      G4double theEkin = boosted.GetKineticEnergy();
      aXsection = ((*G4ParticleHPManager::GetInstance()->GetInelasticFinalStates( dP->GetDefinition() ))[indexEl]
                  ->GetWeightedXsec( theEkin, isotopeJ ) ) * barn;
      if ( aXsection < 0.0 ) {
	if ( failCount < 1000 ) {
	  ++failCount;
          #ifndef G4PHP_DOPPLER_LOOP_ONCE
	  --counter;
	  continue;
          #endif
	} else {
	  aXsection = 0.0;
	}
      }
      // velocity correction
      G4ThreeVector targetVelocity = 1.0 / aThermalNuc.GetMass() * aThermalNuc.GetMomentum();
      aXsection *= (targetVelocity - neutronVelocity).mag() / neutronVMag;
      result += aXsection;
    }
#ifndef G4PHP_DOPPLER_LOOP_ONCE
    size += size;
  }
  result /= counter;
#endif
  return result;
}
