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
//-------------------------------------------------------------------
//
//      Geant4 source file 
//
//      File name: G4ParticleHPIsoProbabilityTable_NJOY.cc
//
//      Authors: Marek Zmeskal (CTU, Czech Technical University in Prague, Czech Republic)
//               Loic Thulliez (CEA France)      
//
//      Creation date: 4 June 2024
//
//      Description: Class for the probability table of the given isotope
//                   and for the given temperature generated with NJOY.
//                   It reads the files with probability tables and
//                   finds the correct cross-section.
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//

#include "G4ParticleHPIsoProbabilityTable_NJOY.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleHPVector.hh"
#include "Randomize.hh"
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
G4ParticleHPIsoProbabilityTable_NJOY::G4ParticleHPIsoProbabilityTable_NJOY() {
  tableOrder = 0;
  lssf_flag = -1;
}

///--------------------------------------------------------------------------------------
G4ParticleHPIsoProbabilityTable_NJOY::~G4ParticleHPIsoProbabilityTable_NJOY() {}

///--------------------------------------------------------------------------------------
void G4ParticleHPIsoProbabilityTable_NJOY::Init( G4int theZ, G4int theA, G4int them, G4double theT, G4String dirName ) {
  Z = theZ;
  A = theA;
  m = them;
  T = theT;  
  G4cout << "The NJOY probability tables are being initialized for Z=" << Z << " A=" << A << " and T=" << T << " K." << G4endl;
  std::string strZ = std::to_string(Z);
  std::string strA = std::to_string(A);
  filename = strZ + "_" + strA;
  if ( m != 0 ) {
    std::string strm = std::to_string(m);
    filename += "_m" + strm;
  }
  G4String fullPathFileName = dirName + filename + "." + std::to_string( (G4int)(T) ) + ".pt.z";
  std::istringstream theData( std::ios::in );
  G4ParticleHPManager::GetInstance()->GetDataStream( fullPathFileName, theData );
  if ( theData.good() ) {    
    G4double emin;
    G4double emax;
    G4double energy;
    theData >> emin >> emax;
    Emin = emin * eV;
    Emax = emax * eV;
    theData >> nEnergies >> tableOrder >> lssf_flag;
    theEnergies = new G4ParticleHPVector(nEnergies);
    theProbabilities = new std::vector< std::vector< G4double >* >;
    theElasticData = new std::vector< std::vector< G4double >* >;
    theCaptureData = new std::vector< std::vector< G4double >* >;
    theFissionData = new std::vector< std::vector< G4double >* >;
    G4double probability, total, elastic, capture, fission;
    for ( G4int i = 0; i < nEnergies; i++ ) {
      theData >> energy;
      theEnergies->SetEnergy( i, energy * eV );
      std::vector< G4double >* vecprob = new std::vector< G4double >;
      std::vector< G4double >* vecela  = new std::vector< G4double >;
      std::vector< G4double >* veccap  = new std::vector< G4double >;
      std::vector< G4double >* vecfiss = new std::vector< G4double >;
      for ( G4int j = 0; j < tableOrder; j++ ) {
        theData >> probability >> total >> elastic >> capture >> fission;
        vecprob->push_back( probability );
        vecela->push_back( elastic );
        veccap->push_back( capture );
        vecfiss->push_back( fission );
      }
      theProbabilities->push_back( vecprob );
      theElasticData->push_back( vecela );
      theCaptureData->push_back( veccap );
      theFissionData->push_back( vecfiss );
    }
    G4cout << "Probability tables found and succesfully read from " << Emin / keV << " keV to " << Emax / keV << " keV." << G4endl;
  } else {
    G4cout << "No probability tables found for this isotope and temperature, smooth cross section will be used instead." << G4endl;
  }
}

///--------------------------------------------------------------------------------------
G4double G4ParticleHPIsoProbabilityTable_NJOY::GetCorrelatedIsoCrossSectionPT( const G4DynamicParticle* dp, G4int MTnumber, 
  const G4Element* ele, G4double& kineticEnergy, G4double& random_number, std::thread::id& id ) {
  if ( kineticEnergy == energy_cache[id] ) {
    if ( MTnumber == 2 ) {           // elastic cross section
      return xsela_cache[id];
    } else if ( MTnumber == 102 ) {  // radiative capture cross section
      return xscap_cache[id];
    } else if ( MTnumber == 18 ) {   // fission cross section
      return xsfiss_cache[id];
    }
  }
  energy_cache[id] = kineticEnergy;
  if ( kineticEnergy < Emin  ||  kineticEnergy > Emax ) {
    // if the kinetic energy is outside of the URR limits for the given isotope, it finds the smooth cross section
    G4int indexEl = (G4int)ele->GetIndex();
    G4int isotopeJ = -1;  // index of isotope in the given element
    for ( G4int j = 0; j < (G4int)ele->GetNumberOfIsotopes(); j++ ) {
      if ( A == (G4int)ele->GetIsotope(j)->GetN() ) {
        isotopeJ = j;
        break;
      }
    }
    G4double frac = ele->GetRelativeAbundanceVector()[isotopeJ];
    G4double weightedelasticXS;
    G4double weightedcaptureXS;
    if ( G4ParticleHPManager::GetInstance()->GetNeglectDoppler() ) {
      weightedelasticXS = (*G4ParticleHPManager::GetInstance()->GetElasticFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
      weightedcaptureXS = (*G4ParticleHPManager::GetInstance()->GetCaptureFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
    } else {
      weightedelasticXS = this->GetDopplerBroadenedElasticXS( dp, indexEl, isotopeJ );
      weightedcaptureXS = this->GetDopplerBroadenedCaptureXS( dp, indexEl, isotopeJ );
    }
    xsela_cache[id] = weightedelasticXS / frac;
    xscap_cache[id] = weightedcaptureXS / frac;
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
  } else if ( lssf_flag == 0 ) {
    G4int indexE = theEnergies->GetEnergyIndex( kineticEnergy );
    std::vector< G4double >* theProbability1 = theProbabilities->at(indexE - 1);
    std::vector< G4double >* theProbability2 = theProbabilities->at(indexE);
    G4double ene1 = theEnergies->GetEnergy(indexE - 1);
    G4double ene2 = theEnergies->GetEnergy(indexE);
    G4double rand = random_number;
    G4int indexP1;
    G4int indexP2;
    for ( indexP1 = 0; indexP1 < tableOrder; indexP1++ ) {
      if ( rand <= theProbability1->at(indexP1) ) break;
    }
    for ( indexP2 = 0; indexP2 < tableOrder; indexP2++ ) {
      if ( rand <= theProbability2->at(indexP2) ) break;
    }
    G4double ela1 = theElasticData->at(indexE - 1)->at(indexP1);
    G4double ela2 = theElasticData->at(indexE)->at(indexP2);
    G4double cap1 = theCaptureData->at(indexE - 1)->at(indexP1);
    G4double cap2 = theCaptureData->at(indexE)->at(indexP2);
    xsela_cache[id] = theInt.Lin( kineticEnergy, ene1, ene2, ela1, ela2 ) * barn;
    xscap_cache[id] = theInt.Lin( kineticEnergy, ene1, ene2, cap1, cap2 ) * barn;
    if ( Z < 88 ) {
      xsfiss_cache[id] = 0.0;
    } else {
      G4double fiss1 = theFissionData->at(indexE - 1)->at(indexP1);
      G4double fiss2 = theFissionData->at(indexE)->at(indexP2);
      xsfiss_cache[id] = theInt.Lin( kineticEnergy, ene1, ene2, fiss1, fiss2 ) * barn;
    }
  } else if ( lssf_flag == 1 ) {
    G4int indexE = theEnergies->GetEnergyIndex( kineticEnergy );
    std::vector< G4double >* theProbability1 = theProbabilities->at(indexE - 1);
    std::vector< G4double >* theProbability2 = theProbabilities->at(indexE);
    G4double ene1 = theEnergies->GetEnergy(indexE - 1);
    G4double ene2 = theEnergies->GetEnergy(indexE);
    G4double rand = random_number;
    G4int indexP1;
    G4int indexP2;
    for ( indexP1 = 0; indexP1 < tableOrder; indexP1++ ) {
      if ( rand <= theProbability1->at(indexP1) ) break;
    }
    for ( indexP2 = 0; indexP2 < tableOrder; indexP2++ ) {
      if ( rand <= theProbability2->at(indexP2) ) break;
    }
    G4int indexEl = (G4int)ele->GetIndex();
    G4int isotopeJ = -1;  // index of isotope in the given element
    for ( G4int j = 0; j < (G4int)ele->GetNumberOfIsotopes(); j++ ) {
      if ( A == (G4int)ele->GetIsotope(j)->GetN() ) {
        isotopeJ = j;
        break;
      }
    }
    G4double frac = ele->GetRelativeAbundanceVector()[isotopeJ];
    G4double weightedelasticXS;
    G4double weightedcaptureXS;
    if ( G4ParticleHPManager::GetInstance()->GetNeglectDoppler() ) {
      weightedelasticXS = (*G4ParticleHPManager::GetInstance()->GetElasticFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
      weightedcaptureXS = (*G4ParticleHPManager::GetInstance()->GetCaptureFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
    } else {
      weightedelasticXS = this->GetDopplerBroadenedElasticXS( dp, indexEl, isotopeJ );
      weightedcaptureXS = this->GetDopplerBroadenedCaptureXS( dp, indexEl, isotopeJ );
    }
    G4double ela1 = theElasticData->at(indexE - 1)->at(indexP1);
    G4double ela2 = theElasticData->at(indexE)->at(indexP2);
    G4double cap1 = theCaptureData->at(indexE - 1)->at(indexP1);
    G4double cap2 = theCaptureData->at(indexE)->at(indexP2);
    G4double elasticXS = weightedelasticXS / frac;
    G4double captureXS = weightedcaptureXS / frac;
    xsela_cache[id] = theInt.Lin( kineticEnergy, ene1, ene2, ela1, ela2 ) * elasticXS;
    xscap_cache[id] = theInt.Lin( kineticEnergy, ene1, ene2, cap1, cap2 ) * captureXS;
    if ( Z < 88 ) {
      xsfiss_cache[id] = 0.0;
    } else {
      if ( G4ParticleHPManager::GetInstance()->GetNeglectDoppler() ) {
        G4double weightedfissionXS = (*G4ParticleHPManager::GetInstance()->GetFissionFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
        G4double fissionXS = weightedfissionXS / frac;
        G4double fiss1 = theFissionData->at(indexE - 1)->at(indexP1);
        G4double fiss2 = theFissionData->at(indexE)->at(indexP2);
        xsfiss_cache[id] = theInt.Lin( kineticEnergy, ene1, ene2, fiss1, fiss2 ) * fissionXS;
      } else {
        G4double weightedfissionXS = this->GetDopplerBroadenedFissionXS( dp, indexEl, isotopeJ );
        G4double fissionXS = weightedfissionXS / frac;
        G4double fiss1 = theFissionData->at(indexE - 1)->at(indexP1);
        G4double fiss2 = theFissionData->at(indexE)->at(indexP2);
        xsfiss_cache[id] = theInt.Lin( kineticEnergy, ene1, ene2, fiss1, fiss2 ) * fissionXS;
      }
    }
  }
  if ( MTnumber == 2 ) {           // elastic cross section
    return xsela_cache[id];
  } else if ( MTnumber == 102 ) {  // radiative capture cross section
    return xscap_cache[id];
  } else if ( MTnumber == 18 ) {   // fission cross section
    return xsfiss_cache[id];
  } else {
    G4cout << "Reaction was not found, returns 0." << G4endl;
    return 0;
  }
}

///--------------------------------------------------------------------------------------
G4double G4ParticleHPIsoProbabilityTable_NJOY::GetIsoCrossSectionPT( const G4DynamicParticle* dp, G4int MTnumber, 
  const G4Element* ele, G4double& kineticEnergy, std::map< std::thread::id, G4double >& random_number_cache, std::thread::id& id ) {
  energy_cache[id] = kineticEnergy;
  if ( kineticEnergy < Emin  ||  kineticEnergy > Emax ) {
    // if the kinetic energy is outside of the URR limits for the given isotope, it finds the smooth cross section
    G4int indexEl = (G4int)ele->GetIndex();
    G4int isotopeJ = -1;  // index of isotope in the given element
    for ( G4int j = 0; j < (G4int)ele->GetNumberOfIsotopes(); j++ ) {
      if ( A == (G4int)ele->GetIsotope(j)->GetN() ) {
        isotopeJ = j;
        break;
      }
    }
    G4double frac = ele->GetRelativeAbundanceVector()[isotopeJ];
    G4double weightedelasticXS;
    G4double weightedcaptureXS;
    if ( G4ParticleHPManager::GetInstance()->GetNeglectDoppler() ) {
      weightedelasticXS = (*G4ParticleHPManager::GetInstance()->GetElasticFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
      weightedcaptureXS = (*G4ParticleHPManager::GetInstance()->GetCaptureFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
    } else {
      weightedelasticXS = this->GetDopplerBroadenedElasticXS( dp, indexEl, isotopeJ );
      weightedcaptureXS = this->GetDopplerBroadenedCaptureXS( dp, indexEl, isotopeJ );
    }
    xsela_cache[id] = weightedelasticXS / frac;
    xscap_cache[id] = weightedcaptureXS / frac;
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
  } else if ( lssf_flag == 0 ) {
    G4int indexE = theEnergies->GetEnergyIndex(kineticEnergy);
    std::vector< G4double >* theProbability1 = theProbabilities->at(indexE - 1);
    std::vector< G4double >* theProbability2 = theProbabilities->at(indexE);
    G4double ene1 = theEnergies->GetEnergy(indexE - 1);
    G4double ene2 = theEnergies->GetEnergy(indexE);
    G4double rand = G4UniformRand();
    random_number_cache[id] = rand;
    G4int indexP1;
    G4int indexP2;
    for ( indexP1 = 0; indexP1 < tableOrder; indexP1++ ) {
      if ( rand <= theProbability1->at(indexP1) ) break;
    }
    for ( indexP2 = 0; indexP2 < tableOrder; indexP2++ ) {
      if ( rand <= theProbability2->at(indexP2) ) break;
    }
    G4double ela1 = theElasticData->at(indexE - 1)->at(indexP1);
    G4double ela2 = theElasticData->at(indexE)->at(indexP2);
    G4double cap1 = theCaptureData->at(indexE - 1)->at(indexP1);
    G4double cap2 = theCaptureData->at(indexE)->at(indexP2);
    xsela_cache[id] = theInt.Lin( kineticEnergy, ene1, ene2, ela1, ela2 ) * barn;
    xscap_cache[id] = theInt.Lin( kineticEnergy, ene1, ene2, cap1, cap2 ) * barn;
    if ( Z < 88 ) {
      xsfiss_cache[id] = 0.0;
    } else {
      G4double fiss1 = theFissionData->at(indexE - 1)->at(indexP1);
      G4double fiss2 = theFissionData->at(indexE)->at(indexP2);
      xsfiss_cache[id] = theInt.Lin( kineticEnergy, ene1, ene2, fiss1, fiss2 ) * barn;
    }
  } else if ( lssf_flag == 1 ) {
    G4int indexE = theEnergies->GetEnergyIndex( kineticEnergy );
    std::vector< G4double >* theProbability1 = theProbabilities->at(indexE - 1);
    std::vector< G4double >* theProbability2 = theProbabilities->at(indexE);
    G4double ene1 = theEnergies->GetEnergy(indexE - 1);
    G4double ene2 = theEnergies->GetEnergy(indexE);
    G4double rand = G4UniformRand();
    random_number_cache[id] = rand;
    G4int indexP1;
    G4int indexP2;
    for ( indexP1 = 0; indexP1 < tableOrder; indexP1++ ) {
      if ( rand <= theProbability1->at(indexP1) ) break;
    }
    for ( indexP2 = 0; indexP2 < tableOrder; indexP2++ ) {
      if ( rand <= theProbability2->at(indexP2) ) break;
    }
    G4int indexEl = (G4int)ele->GetIndex();
    G4int isotopeJ = -1;  // index of isotope in the given element
    for ( G4int j = 0; j < (G4int)ele->GetNumberOfIsotopes(); j++ ) {
      if ( A == (G4int)ele->GetIsotope(j)->GetN() ) {
        isotopeJ = j;
        break;
      }
    }
    G4double frac = ele->GetRelativeAbundanceVector()[isotopeJ];
    G4double weightedelasticXS;
    G4double weightedcaptureXS;
    if ( G4ParticleHPManager::GetInstance()->GetNeglectDoppler() ) {
      weightedelasticXS = (*G4ParticleHPManager::GetInstance()->GetElasticFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
      weightedcaptureXS = (*G4ParticleHPManager::GetInstance()->GetCaptureFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
    } else {
      weightedelasticXS = this->GetDopplerBroadenedElasticXS( dp, indexEl, isotopeJ );
      weightedcaptureXS = this->GetDopplerBroadenedCaptureXS( dp, indexEl, isotopeJ );
    }
    G4double ela1 = theElasticData->at(indexE - 1)->at(indexP1);
    G4double ela2 = theElasticData->at(indexE)->at(indexP2);
    G4double cap1 = theCaptureData->at(indexE - 1)->at(indexP1);
    G4double cap2 = theCaptureData->at(indexE)->at(indexP2);
    G4double elasticXS = weightedelasticXS / frac;
    G4double captureXS = weightedcaptureXS / frac;
    xsela_cache[id] = theInt.Lin( kineticEnergy, ene1, ene2, ela1, ela2 ) * elasticXS;
    xscap_cache[id] = theInt.Lin( kineticEnergy, ene1, ene2, cap1, cap2 ) * captureXS;
    if ( Z < 88 ) {
      xsfiss_cache[id] = 0.0;
    } else {
      if ( G4ParticleHPManager::GetInstance()->GetNeglectDoppler() ) {
          G4double weightedfissionXS = (*G4ParticleHPManager::GetInstance()->GetFissionFinalStates())[indexEl]->GetWeightedXsec( kineticEnergy, isotopeJ );
          G4double fissionXS = weightedfissionXS / frac;
          G4double fiss1 = theFissionData->at(indexE - 1)->at(indexP1);
          G4double fiss2 = theFissionData->at(indexE)->at(indexP2);
          xsfiss_cache[id] = theInt.Lin(kineticEnergy, ene1, ene2, fiss1, fiss2) * fissionXS;
      } else {
          G4double weightedfissionXS = this->GetDopplerBroadenedFissionXS( dp, indexEl, isotopeJ );
          G4double fissionXS = weightedfissionXS / frac;
          G4double fiss1 = theFissionData->at(indexE - 1)->at(indexP1);
          G4double fiss2 = theFissionData->at(indexE)->at(indexP2);
          xsfiss_cache[id] = theInt.Lin( kineticEnergy, ene1, ene2, fiss1, fiss2 ) * fissionXS;
      }
    }
  }
  if ( MTnumber == 2 ) {           // elastic cross section
    return xsela_cache[id];
  } else if ( MTnumber == 102 ) {  // radiative capture cross section
    return xscap_cache[id];
  } else if ( MTnumber == 18 ) {   // fission cross section
    return xsfiss_cache[id];
  } else {
    G4cout << "Reaction was not found, returns 0." << G4endl;
    return 0;
  }
}
