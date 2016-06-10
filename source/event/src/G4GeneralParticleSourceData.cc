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
// GEANT4 Class file
//
//
// File name:     G4GeneralParticleSourceData.cc
//
// Author:        Andrew Green
//
// Creation date: 20 Mar 2014
//
// Modifications:
// 24/03/2014
// Fixed a bug whereby there was a data race for ownership of the "currentSource"
// member data. This has been resolved by returning a pointer to the G4SPS from
// the source vector. Tested up to 4 threads, and works fine; may need some further
// locking in the G4SPS code.
//
//
// Class Description:
// This class uses the singleton pattern to create a single copy of the data
// needed for the G4GPS class. As yet, only the largest parts have been split
// off.
//
//

#include "G4GeneralParticleSourceData.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"




namespace
{
  G4Mutex singMutex = G4MUTEX_INITIALIZER; //Protects singleton access
}

//G4GeneralParticleSourceData* G4GeneralParticleSourceData::theInstance = 0;

G4GeneralParticleSourceData::G4GeneralParticleSourceData() :
		multiple_vertex(false) ,flat_sampling(false),
		normalised(false),currentSourceIdx(0)
{
    G4MUTEXINIT(mutex);
    
    sourceVector.clear();
    sourceIntensity.clear();
    sourceProbability.clear();
    
    currentSource = new G4SingleParticleSource();
    sourceVector.push_back(currentSource);
    sourceIntensity.push_back(1.);

}

G4GeneralParticleSourceData::~G4GeneralParticleSourceData()
{
  G4MUTEXDESTROY(mutex);
  for ( std::vector<G4SingleParticleSource*>::const_iterator it = sourceVector.begin() ;
	it != sourceVector.end() ; ++it ) { delete *it; }
  sourceVector.clear();
}


G4GeneralParticleSourceData* G4GeneralParticleSourceData::Instance()
{
    G4AutoLock lock(&singMutex);
    static G4GeneralParticleSourceData instance;
    return &instance;
}

void G4GeneralParticleSourceData::IntensityNormalise()
{
    G4double total  = 0.;
    size_t i = 0 ;
    for (i = 0; i < sourceIntensity.size(); i++)
    {
        total += sourceIntensity[i] ;
    }
    sourceProbability.clear();
    std::vector <G4double> sourceNormalizedIntensity;
    sourceNormalizedIntensity.clear();
    
    sourceNormalizedIntensity.push_back(sourceIntensity[0]/total);
    sourceProbability.push_back(sourceNormalizedIntensity[0]);

    for ( i = 1 ;  i < sourceIntensity.size(); i++)
    {
        sourceNormalizedIntensity.push_back(sourceIntensity[i]/total);
        sourceProbability.push_back(sourceNormalizedIntensity[i] + sourceProbability[i-1]);
    }

    // set source weights here based on sampling scheme (analog/flat) and intensities
    for ( i = 0 ;  i < sourceIntensity.size(); i++)
    {
        if (!flat_sampling)
        {
            this->GetCurrentSource(i)->GetBiasRndm()->SetIntensityWeight(1.);
        }
        else
        {
            this->GetCurrentSource(i)->GetBiasRndm()->SetIntensityWeight(sourceNormalizedIntensity[i]*sourceIntensity.size());
        }
    }

    normalised = true;
}

void G4GeneralParticleSourceData::SetCurrentSourceIntensity(G4double intensity)
{
    sourceIntensity.at(currentSourceIdx) = intensity;
    normalised = false;
}

void G4GeneralParticleSourceData::AddASource(G4double intensity)
{
    currentSource = new G4SingleParticleSource();
    sourceVector.push_back(currentSource);
    sourceIntensity.push_back(intensity);
    currentSourceIdx = sourceVector.size() - 1;
    normalised = false;
}

void G4GeneralParticleSourceData::DeleteASource(G4int idx)
{
    delete sourceVector[idx];
    sourceVector.erase(sourceVector.begin() + idx);
    sourceIntensity.erase(sourceIntensity.begin()+idx);
    normalised = false ;
    if (currentSourceIdx == idx )
    {
        if ( this->GetIntensityVectorSize() > 0 )
        {
            currentSource = this->GetCurrentSource(0);
            currentSourceIdx = 0;
        }
        else
        {
            currentSource = NULL;
            currentSourceIdx = -1;
        }
    }

}

void G4GeneralParticleSourceData::ClearSources()
{
    currentSourceIdx = -1;
    currentSource = NULL;
    for ( std::vector<G4SingleParticleSource*>::iterator it = sourceVector.begin();
    	  it != sourceVector.end() ; ++it ) { delete *it; }
    sourceVector.clear();
    sourceIntensity.clear();
    normalised = false;
}

void G4GeneralParticleSourceData::SetVerbosityAllSources(G4int vl )
{
    for ( std::vector<G4SingleParticleSource*>::iterator it = sourceVector.begin();
    	  it != sourceVector.end() ; ++it ) {
    	(*it)->SetVerbosity(vl);

    }

}

G4SingleParticleSource* G4GeneralParticleSourceData::GetCurrentSource(G4int idx)
{
    currentSource = sourceVector[idx];
    currentSourceIdx = idx;
    return sourceVector[idx];
}

void G4GeneralParticleSourceData::Lock()
{
    G4MUTEXLOCK(&mutex);
}

void G4GeneralParticleSourceData::Unlock()
{
    G4MUTEXUNLOCK(&mutex);
}
