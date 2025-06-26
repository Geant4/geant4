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
//      File name: G4ParticleHPIsoProbabilityTable_CALENDF.cc
//
//      Authors: Marek Zmeskal (CTU, Czech Technical University in Prague, Czech Republic)
//               Loic Thulliez (CEA France)      
//
//      Creation date: 4 June 2024
//
//      Description: Class for the probability table of the given isotope
//                   and for the given temperature generated with CALENDF.
//                   It reads the files with probability tables and
//                   finds the correct cross-section.
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//

#include "G4ParticleHPIsoProbabilityTable_CALENDF.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ParticleHPVector.hh"
#include "G4DynamicParticle.hh"
#include "G4NucleiProperties.hh"
#include "G4ReactionProduct.hh"
#include "G4ParticleHPManager.hh"
#include "G4ParticleHPChannel.hh"
#include "G4ParticleHPChannelList.hh"
#include "G4Nucleus.hh"
#include "G4Element.hh"
#include <string>
#include <sstream>

///--------------------------------------------------------------------------------------
G4ParticleHPIsoProbabilityTable_CALENDF::G4ParticleHPIsoProbabilityTable_CALENDF() : theInelasticData(nullptr) {}

///--------------------------------------------------------------------------------------
G4ParticleHPIsoProbabilityTable_CALENDF::~G4ParticleHPIsoProbabilityTable_CALENDF() {
  for ( auto it = theInelasticData->cbegin(); it != theInelasticData->cend(); ++it ) {
    delete* it;
  }
  delete theInelasticData;
}

///--------------------------------------------------------------------------------------
void G4ParticleHPIsoProbabilityTable_CALENDF::Init( G4int theZ, G4int theA, G4int them, G4double theT, const G4String& dirName ) {
  Z = theZ;
  A = theA;
  m = them;
  T = theT;  
  G4cout << "The CALENDF probability tables are being initialized for Z=" << Z << " A=" << A << " and T=" << T << " K." << G4endl;
  filename = std::to_string(Z) + "_" + std::to_string(A);
  if ( m != 0 ) filename += "_m" + std::to_string(m);
  G4String fullPathFileName = dirName + filename + "." + std::to_string( (G4int)(T) ) + ".pt"; 
  std::istringstream theData( std::ios::in );
  G4ParticleHPManager::GetInstance()->GetDataStream( fullPathFileName, theData );
  if ( theData.good() ) {    
    G4double emin;
    G4double emax;
    theData >> emin >> emax;
    Emin = emin * eV;
    Emax = emax * eV;
    theData >> nEnergies;
    theEnergies = new G4ParticleHPVector( nEnergies );
    theProbabilities = new std::vector< std::vector< G4double >* >;
    theElasticData = new std::vector< std::vector< G4double >* >;
    theCaptureData = new std::vector< std::vector< G4double >* >;
    theFissionData = new std::vector< std::vector< G4double >* >;
    theInelasticData = new std::vector< std::vector< G4double >* >;
    G4double tableOrder;
    G4double probability, total, elastic, capture, fission, inelastic;
    for ( G4int i = 0; i < nEnergies; i++ ) {
      theData >> emin >> emax >> tableOrder;
      theEnergies->SetData( i, emax * eV, tableOrder );
      std::vector< G4double >* vecprob = new std::vector< G4double >;
      std::vector< G4double >* vecela = new std::vector< G4double >;
      std::vector< G4double >* veccap = new std::vector< G4double >;
      std::vector< G4double >* vecfiss = new std::vector< G4double >;
      std::vector< G4double >* vecinela = new std::vector< G4double >;
      for ( G4int j = 0; j < tableOrder; j++ ) {
  	theData >> probability >> total >> elastic >> capture >> fission >> inelastic;
  	vecprob->push_back( probability );
  	vecela->push_back( elastic * barn );
  	veccap->push_back( capture * barn );
  	vecfiss->push_back( fission * barn );
  	vecinela->push_back( inelastic * barn );
      }
      theProbabilities->push_back( vecprob );
      theElasticData->push_back( vecela );
      theCaptureData->push_back( veccap );
      theFissionData->push_back( vecfiss );
      theInelasticData->push_back( vecinela );
    }  
    G4cout << "Probability tables found and succesfully read from " << Emin / keV << " keV to " << Emax / keV << " keV." << G4endl;
  } else {
    G4cout << "No probability tables found for this isotope and temperature, smooth cross section will be used instead." << G4endl;
  }
}

///--------------------------------------------------------------------------------------
G4double G4ParticleHPIsoProbabilityTable_CALENDF::GetCorrelatedIsoCrossSectionPT( const G4DynamicParticle* dp, 
  G4int MTnumber, const G4Element* ele, G4double &kineticEnergy, G4double &random_number, std::thread::id &id ) {
  // if this function is called again for different reaction or the particle came to different region 
  // with the same temperature and unchanged temperature
  if ( kineticEnergy == energy_cache[id] ) {
    if ( MTnumber == 2 ) {  // elastic cross section
      return xsela_cache[id];
    } else if ( MTnumber == 102 ) {  // radiative capture cross section
      return xscap_cache[id];
    } else if ( MTnumber == 18 ) {  // fission cross section
      return xsfiss_cache[id];
    } else if ( MTnumber == 3 ) {  // inelastic cross section
      return xsinela_cache[id];
    }
  }
  energy_cache[id] = kineticEnergy;
  if ( kineticEnergy < Emin  ||  kineticEnergy > Emax ) {
    // if the kinetic energy is outside of the URR limits for the given isotope, it finds the smooth cross section
    G4int indexEl = (G4int)ele->GetIndex();
    G4int isotopeJ = 0;  // index of isotope in the given element
    G4int n_isotopes = (G4int)ele->GetNumberOfIsotopes();
    for ( G4int j = 0; j < n_isotopes; j++ ) {
      if ( A == (G4int)( ele->GetIsotope(j)->GetN() ) ) {
	isotopeJ = j;
	break;
      }
    }
    G4double frac = ele->GetRelativeAbundanceVector()[isotopeJ];
    G4double weightedelasticXS;
    G4double weightedcaptureXS;
    G4double weightedinelasticXS;
    if ( G4ParticleHPManager::GetInstance()->GetNeglectDoppler() ) {
      weightedelasticXS = (*G4ParticleHPManager::GetInstance()->GetElasticFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
      weightedcaptureXS = (*G4ParticleHPManager::GetInstance()->GetCaptureFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
      weightedinelasticXS = ((*G4ParticleHPManager::GetInstance()->GetInelasticFinalStates(dp->GetDefinition()))[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ )) * barn;
    } else {
      weightedelasticXS = this->GetDopplerBroadenedElasticXS( dp, indexEl, isotopeJ );
      weightedcaptureXS = this->GetDopplerBroadenedCaptureXS( dp, indexEl, isotopeJ );
      weightedinelasticXS = this->GetDopplerBroadenedInelasticXS( dp, indexEl, isotopeJ );
    }
    xsela_cache[id] = weightedelasticXS / frac;
    xscap_cache[id] = weightedcaptureXS / frac;
    xsinela_cache[id] = weightedinelasticXS / frac;
    if ( Z < 88 ) {
      xsfiss_cache[id] = 0.0;
    } else {
      if ( G4ParticleHPManager::GetInstance()->GetNeglectDoppler() ) {
	G4double weightedfissionXS = (*G4ParticleHPManager::GetInstance()->GetFissionFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
	xsfiss_cache[id] = weightedfissionXS / frac;
      } else {
	G4double weightedfissionXS = this->GetDopplerBroadenedFissionXS( dp, indexEl, isotopeJ );
	xsfiss_cache[id] = weightedfissionXS / frac;
      }
    }
  } else {
    G4int indexE = theEnergies->GetEnergyIndex( kineticEnergy );
    G4int order = (G4int)( theEnergies->GetY(indexE) + 0.5 );
    std::vector< G4double >* theProbability = theProbabilities->at( indexE );
    G4double rand = random_number;
    G4int indexP;
    for ( indexP = 0; indexP < order; indexP++ ) {
      if ( rand <= theProbability->at(indexP) ) break;
    }
    xsela_cache[id] = theElasticData->at(indexE)->at(indexP);
    xscap_cache[id] = theCaptureData->at(indexE)->at(indexP);
    xsinela_cache[id] = theInelasticData->at(indexE)->at(indexP);
    if ( Z < 88 ) xsfiss_cache[id] = 0.0;
    else          xsfiss_cache[id] = theFissionData->at(indexE)->at(indexP);
  }
  if ( MTnumber == 2 ) {           // elastic cross section
    return xsela_cache[id];
  } else if ( MTnumber == 102 ) {  // radiative capture cross section
    return xscap_cache[id];
  } else if ( MTnumber == 18 ) {   // fission cross section
    return xsfiss_cache[id];
  } else if ( MTnumber == 3 ) {    // inelastic cross section
    return xsinela_cache[id];
  } else {
    G4cout << "Reaction was not found, returns 0." << G4endl;
    return 0;
  }
}

///--------------------------------------------------------------------------------------
G4double G4ParticleHPIsoProbabilityTable_CALENDF::GetIsoCrossSectionPT( const G4DynamicParticle* dp, G4int MTnumber, 
									const G4Element* ele,G4double &kineticEnergy,
									std::map< std::thread::id, G4double >& random_number_cache,
									std::thread::id &id ) {
  energy_cache[id] = kineticEnergy;
  if ( kineticEnergy < Emin  ||  kineticEnergy > Emax ) {
    // if the kinetic energy is outside of the URR limits for the given isotope, it finds the smooth cross section
    G4int indexEl = (G4int)ele->GetIndex();
    G4int isotopeJ = -1;  //index of isotope in the given element
    G4int n_isotopes = (G4int)ele->GetNumberOfIsotopes();
    for ( G4int j = 0; j < n_isotopes; j++ ) {
      if ( A == (G4int)ele->GetIsotope(j)->GetN() ) {
  	isotopeJ = j;
  	break;
      }
    }
    if (isotopeJ == -1) { return 0.0; }
    G4double frac = ele->GetRelativeAbundanceVector()[isotopeJ]; 
    G4double weightedelasticXS;
    G4double weightedcaptureXS;
    G4double weightedinelasticXS;
    auto manHP = G4ParticleHPManager::GetInstance(); 
    if (manHP->GetNeglectDoppler()) {
      weightedelasticXS = (*manHP->GetElasticFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
      weightedcaptureXS = (*manHP->GetCaptureFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
      weightedinelasticXS = (*manHP->GetInelasticFinalStates(dp->GetDefinition()))[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ ) * barn;
    } else {
      weightedelasticXS = GetDopplerBroadenedElasticXS( dp, indexEl, isotopeJ );
      weightedcaptureXS = GetDopplerBroadenedCaptureXS( dp, indexEl, isotopeJ );
      weightedinelasticXS = GetDopplerBroadenedInelasticXS( dp, indexEl, isotopeJ );
    }
    xsela_cache[id] = weightedelasticXS / frac;
    xscap_cache[id] = weightedcaptureXS / frac;
    xsinela_cache[id] = weightedinelasticXS / frac;
    if ( Z < 88 ) {
      xsfiss_cache[id] = 0.0;
    } else {
      if ( manHP->GetNeglectDoppler() ) {
  	G4double weightedfissionXS = (*manHP->GetFissionFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
  	xsfiss_cache[id] = weightedfissionXS / frac;
      } else {
  	G4double weightedfissionXS = GetDopplerBroadenedFissionXS( dp, indexEl, isotopeJ );
  	xsfiss_cache[id] = weightedfissionXS / frac;
      }
    }
  } else {
    G4int indexE = theEnergies->GetEnergyIndex( kineticEnergy );
    G4int order = (G4int)( theEnergies->GetY(indexE) + 0.5 );
    std::vector< G4double >* theProbability = theProbabilities->at(indexE);
    G4double rand = G4UniformRand();
    random_number_cache[id] = rand;
    G4int indexP;
    for ( indexP = 0; indexP < order; indexP++ ) {
      if ( rand <= theProbability->at(indexP) ) break;
    }
    xsela_cache[id] = theElasticData->at(indexE)->at(indexP);
    xscap_cache[id] = theCaptureData->at(indexE)->at(indexP);
    xsinela_cache[id] = theInelasticData->at(indexE)->at(indexP);
    if ( Z < 88 ) xsfiss_cache[id] = 0.0;
    else          xsfiss_cache[id] = theFissionData->at(indexE)->at(indexP);
  }
  if (MTnumber == 2) {           // elastic cross section
    return xsela_cache[id];
  } else if (MTnumber == 102) {  // radiative capture cross section
    return xscap_cache[id];
  } else if (MTnumber == 18) {   // fission cross section
    return xsfiss_cache[id];
  } else if (MTnumber == 3) {    // inelastic cross section
    return xsinela_cache[id];
  } else {
    //G4cout << "Reaction was not found, returns 0." << G4endl;
    return 0;
  }
}
