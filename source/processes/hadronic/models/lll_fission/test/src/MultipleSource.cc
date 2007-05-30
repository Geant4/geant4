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
// MODULE:       MultipleSource.cc
//
// Version:      1.0
// Date:         1/12/07
// Author:       Jerome M Verbeke 
// Organisation: LLNL
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
// Version 1.0, 01/12/2007, Jerome M Verbeke, Created.
//
//
///////////////////////////////////////////////////////////////////////////////
//
#include "G4Event.hh"
#include "Randomize.hh"
#include "MultipleSource.hh"

MultipleSource::MultipleSource()
  : multiple_vertex(false), flat_sampling(false), weight_change(1.)
{
  sourceVector.clear();
  sourceIntensity.clear();
  sourceProbability.clear();
  currentSource = new SingleSource();
  sourceVector.push_back(currentSource);
  sourceIntensity.push_back(1.);
  currentSourceIdx = G4int(sourceVector.size() - 1);
  IntensityNormalization();
}

MultipleSource::MultipleSource(SingleSource* src, G4double strength)
  : multiple_vertex(false), flat_sampling(false), weight_change(1.)
{
  sourceVector.clear();
  sourceIntensity.clear();
  sourceProbability.clear();
  currentSource = src;
  sourceVector.push_back(currentSource);
  sourceIntensity.push_back(strength);
  currentSourceIdx = G4int(sourceVector.size() - 1);
  normalised = false;
}
  
MultipleSource::~MultipleSource()
{
}

void MultipleSource::AddaSource(G4double aV)
{
  currentSource = new SingleSource();
  sourceVector.push_back(currentSource);
  sourceIntensity.push_back(aV);
  currentSourceIdx = G4int(sourceVector.size() - 1);
  IntensityNormalization();
}

void MultipleSource::AddaSource(SingleSource* src, G4double aV)
{
  currentSource = src;
  sourceVector.push_back(currentSource);
  sourceIntensity.push_back(aV);
  currentSourceIdx = G4int(sourceVector.size() - 1);
  normalised = false;
}

void MultipleSource::IntensityNormalization()
{
  G4double total  = 0.;
  size_t i = 0 ;
  for (i = 0; i < sourceIntensity.size(); i++) 
    total += sourceIntensity[i] ;
  //
  sourceProbability.clear();
  sourceIntensity[0] =  sourceIntensity[0]/total;
  sourceProbability.push_back(sourceIntensity[0]);
  for ( i = 1 ;  i < sourceIntensity.size(); i++) {
    sourceIntensity[i] = sourceIntensity[i]/total;
    sourceProbability.push_back(sourceIntensity[i] + sourceProbability[i-1]);
  } 

  // set source weights here based on sampling scheme (analog/flat) and intensities
  for ( i = 0 ;  i < sourceIntensity.size(); i++) {
    if (!flat_sampling) {
      sourceVector[i]->GetBiasRndm()->SetIntensityWeight(1.);
    } else {
      sourceVector[i]->GetBiasRndm()->SetIntensityWeight(sourceIntensity[i]*sourceIntensity.size());
    }
  }

  normalised = true;
} 

void MultipleSource::ListSource()
{
  G4cout << " The number of particle source is " << sourceIntensity.size() << G4endl;
  for (size_t i = 0 ; i < sourceIntensity.size(); i++) 
    G4cout << "   source " << i << " intensity is " << sourceIntensity[i] << G4endl;
}

void MultipleSource::SetCurrentSourceto(G4int aV)
{
  size_t id = size_t (aV) ;
  if ( id <= sourceIntensity.size() ) {
    currentSourceIdx = aV;
    currentSource = sourceVector[id];
  } else {
    G4cout << " source index is invalid " << G4endl;
    G4cout << "    it shall be <= " << sourceIntensity.size() << G4endl;
  }
}

void MultipleSource::SetCurrentSourceIntensity(G4double aV)
{
  sourceIntensity[currentSourceIdx] = aV;
  normalised = false;
}

void MultipleSource::ClearAll()
{
  currentSourceIdx = -1;
  currentSource = 0;
  sourceVector.clear();
  sourceIntensity.clear();
  sourceProbability.clear();
}

void MultipleSource::DeleteaSource(G4int aV)
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

void MultipleSource::GeneratePrimaryVertex(G4Event* evt)
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
