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
///////////////////////////////////////////////////////////////////////////////
//
// MODULE:       G4GeneralParticleSource.cc
//
// Version:      2.0
// Date:         5/02/04
// Author:       Fan Lei 
// Organisation: QinetiQ ltd.
// Customer:     ESA/ESTEC
//
// Documentation avaialable at http://reat.space.qinetiq.com/gps
//   These include:
//       User Requirement Document (URD)
//       Software Specification Documents (SSD)
//       Software User Manual (SUM): on-line version available
//       Technical Note (TN) on the physics and algorithms
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
// Version 2.0, 05/02/2004, Fan Lei, Created.
//    After changes to version 1.1 as in Geant4 v6.0
//     - Mutilple particle source definition
//     - Re-structured commands
//     - Split the task into smaller classes
//
//     - old commonds have been retained for backward compatibility, will be
//       removed in the future. 
//
//
///////////////////////////////////////////////////////////////////////////////
//
#include "G4Event.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"

G4GeneralParticleSource::G4GeneralParticleSource()
  : multiple_vertex(false), flat_sampling(false), weight_change(1.)
{
  sourceVector.clear();
  sourceIntensity.clear();
  sourceProbability.clear();
  currentSource = new G4SingleParticleSource();
  sourceVector.push_back(currentSource);
  sourceIntensity.push_back(1.);
  currentSourceIdx = G4int(sourceVector.size() - 1);
  theMessenger = new G4GeneralParticleSourceMessenger(this);
  theMessenger->SetParticleGun(currentSource);
  IntensityNormalization();
}
  
G4GeneralParticleSource::~G4GeneralParticleSource()
{
  delete theMessenger;
}

void G4GeneralParticleSource::AddaSource(G4double aV)
{
  currentSource = new G4SingleParticleSource();
  theMessenger->SetParticleGun(currentSource);
  sourceVector.push_back(currentSource);
  sourceIntensity.push_back(aV);
  currentSourceIdx = G4int(sourceVector.size() - 1);
  IntensityNormalization();
}

void G4GeneralParticleSource::IntensityNormalization()
{
  G4double total  = 0.;
  size_t i = 0 ;
  for (i = 0; i < sourceIntensity.size(); i++) 
    total += sourceIntensity[i] ;
  //
  sourceProbability.clear();
  std::vector <G4double> sourceNormalizedIntensity;
  sourceNormalizedIntensity.clear();

  sourceNormalizedIntensity.push_back(sourceIntensity[0]/total);
  sourceProbability.push_back(sourceNormalizedIntensity[0]);

  for ( i = 1 ;  i < sourceIntensity.size(); i++) {
    sourceNormalizedIntensity.push_back(sourceIntensity[i]/total);
    sourceProbability.push_back(sourceNormalizedIntensity[i] + sourceProbability[i-1]);
  }

  // set source weights here based on sampling scheme (analog/flat) and intensities
  for ( i = 0 ;  i < sourceIntensity.size(); i++) {
    if (!flat_sampling) {
      sourceVector[i]->GetBiasRndm()->SetIntensityWeight(1.);
    } else {
      sourceVector[i]->GetBiasRndm()->SetIntensityWeight(sourceNormalizedIntensity[i]*sourceIntensity.size());
    }
  }

  normalised = true;
} 

void G4GeneralParticleSource::ListSource()
{
  G4cout << " The number of particle sources is " << sourceIntensity.size() << G4endl;
  for (size_t i = 0 ; i < sourceIntensity.size(); i++)
    G4cout << "   source " << i << " intensity is " << sourceIntensity[i] << G4endl;
}

void G4GeneralParticleSource::SetCurrentSourceto(G4int aV)
{
  size_t id = size_t (aV) ;
  if ( id <= sourceIntensity.size() ) {
    currentSourceIdx = aV;
    currentSource = sourceVector[id];
    theMessenger->SetParticleGun(currentSource);
    //
  } else {
    G4cout << " source index is invalid " << G4endl;
    G4cout << "    it shall be <= " << sourceIntensity.size() << G4endl;
  }
}

void G4GeneralParticleSource::SetCurrentSourceIntensity(G4double aV)
{
  sourceIntensity[currentSourceIdx] = aV;
  normalised = false;
}

void G4GeneralParticleSource::ClearAll()
{
  currentSourceIdx = -1;
  currentSource = 0;
  sourceVector.clear();
  sourceIntensity.clear();
  sourceProbability.clear();
}

void G4GeneralParticleSource::DeleteaSource(G4int aV)
{
  size_t id = size_t (aV) ;
  if ( id <= sourceIntensity.size() ) {
    sourceVector.erase(sourceVector.begin()+aV);
    sourceIntensity.erase(sourceIntensity.begin()+aV);
    normalised = false ;
    if (currentSourceIdx == aV ) { 
	if ( sourceIntensity.size() > 0 ) { 
	  currentSource = sourceVector[0];
	  currentSourceIdx = 1;
	} else {
	  currentSource = 0;
	  currentSourceIdx = -1;
	}
    }	  		
  } else {
    G4cout << " source index is invalid " << G4endl;
    G4cout << "    it shall be <= " << sourceIntensity.size() << G4endl;
  }
} 

void G4GeneralParticleSource::GeneratePrimaryVertex(G4Event* evt)
{
  if (!multiple_vertex){
    if (sourceIntensity.size() > 1) {
      if (!normalised) IntensityNormalization();
      G4double rndm = G4UniformRand();
      size_t i = 0 ;
      if (!flat_sampling) {
	while ( rndm > sourceProbability[i] ) i++;
	(currentSource = sourceVector[i]);
      } else {
	i = size_t (sourceIntensity.size()*rndm);
	currentSource = sourceVector[i];
      }
    }
    currentSource-> GeneratePrimaryVertex(evt);
  } 
  else {
    for (size_t i = 0; i <  sourceIntensity.size(); i++) {
      sourceVector[i]->GeneratePrimaryVertex(evt); 
    }
  }
}
