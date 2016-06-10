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
//  25/03/2014, Andrew Green
//      Various changes to use the new G4GeneralParticleSourceData class, mostly
//      just transparent wrappers around the thread safe object.
//
///////////////////////////////////////////////////////////////////////////////
//
#include "G4Event.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"
#include "G4SingleParticleSource.hh"

#include "G4GeneralParticleSourceData.hh"

#include "G4Threading.hh"
#include "G4AutoLock.hh"

namespace {
    G4Mutex messangerInit = G4MUTEX_INITIALIZER;
}

G4GeneralParticleSource::G4GeneralParticleSource() : multiple_vertex(false), flat_sampling(false),normalised(false),
    theMessenger(0)
{
    GPSData = G4GeneralParticleSourceData::Instance();
    currentSource = GPSData->GetCurrentSource();
    currentSourceIdx = G4int(GPSData->GetSourceVectorSize() - 1);

    //Messenger is special, only a worker should instantiate it. Singleton pattern
    theMessenger = G4GeneralParticleSourceMessenger::GetInstance(this);
    //Some initialization should be done only once
    G4AutoLock l(&messangerInit);
    static G4bool onlyOnce = false;
    if ( !onlyOnce ) {
        theMessenger->SetParticleGun(currentSource);
        IntensityNormalization();
        onlyOnce = true;
    }

}

G4GeneralParticleSource::~G4GeneralParticleSource()
{
    theMessenger->Destroy();
}

void G4GeneralParticleSource::AddaSource(G4double aV)
{
    normalised=false;
    GPSData->Lock();
    GPSData->AddASource(aV);
    currentSource = GPSData->GetCurrentSource();
    theMessenger->SetParticleGun(currentSource);
    currentSourceIdx = G4int(GPSData->GetSourceVectorSize() - 1);
    IntensityNormalization();
    GPSData->Unlock();
}

void G4GeneralParticleSource::IntensityNormalization()
{
    GPSData->IntensityNormalise();
    normalised=true;
}

void G4GeneralParticleSource::ListSource()
{
    G4cout << "The number of particle sources is: " << GPSData->GetIntensityVectorSize() << G4endl;
    for(G4int i=0; i<GPSData->GetIntensityVectorSize(); i++)
    {
        G4cout << "\tsource " << i << " intensity is: " << GPSData->GetIntensity(i) << G4endl;;
    }
}

void G4GeneralParticleSource::SetCurrentSourceto(G4int aV)
{
    G4int id = aV;
    if ( id <= GPSData->GetIntensityVectorSize() )
    {
        currentSourceIdx = aV;
        currentSource = GPSData->GetCurrentSource(id);
        theMessenger->SetParticleGun(currentSource);
    }
    else
    {
        G4cout << " source index is invalid " << G4endl;
        G4cout << "    it shall be <= " << GPSData->GetIntensityVectorSize() << G4endl;
    }
}

void G4GeneralParticleSource::SetCurrentSourceIntensity(G4double aV)
{
    GPSData->Lock();
    GPSData->SetCurrentSourceIntensity(aV);
    GPSData->Unlock();
    normalised = false;
}

void G4GeneralParticleSource::ClearAll()
{
    currentSourceIdx = -1;
    currentSource = 0;
    GPSData->ClearSources();
    normalised=false;
}

void G4GeneralParticleSource::DeleteaSource(G4int aV)
{
    G4int id = aV;
    if ( id <= GPSData->GetIntensityVectorSize() )
    {
        GPSData->DeleteASource(aV);
        normalised=false;
    }
    else
    {
        G4cout << " source index is invalid " << G4endl;
        G4cout << "    it shall be <= " << GPSData->GetIntensityVectorSize() << G4endl;
    }
}

void G4GeneralParticleSource::GeneratePrimaryVertex(G4Event* evt)
{
    if (!multiple_vertex)
    {
        if (GPSData->GetIntensityVectorSize() > 1)
        {
            //Try to minimize locks
            if (! normalised ) {
                //According to local variable, normalization is needed
                //Check with underlying (shared resource), another
                //thread could have already normalized this
                GPSData->Lock();
                G4bool norm = GPSData->Normalised();
                if (!norm) {
                    IntensityNormalization();
                }
                //This takes care of the case in which the local variable
                //is False and the underlying source is.
                normalised = GPSData->Normalised();
                GPSData->Unlock();
            }
            G4double rndm = G4UniformRand();
            size_t i = 0 ;
            if (!flat_sampling)
            {
                while ( rndm > GPSData->GetSourceProbability(i) ) i++;
                    (currentSource = GPSData->GetCurrentSource(i));
            }
            else
            {
                i = size_t (GPSData->GetIntensityVectorSize()*rndm);
                currentSource = GPSData->GetCurrentSource(i);
            }
        }
        currentSource-> GeneratePrimaryVertex(evt);
    } 
    else
    {
        for (G4int i = 0; i <  GPSData->GetIntensityVectorSize(); i++)
        {
            GPSData->GetCurrentSource(i)->GeneratePrimaryVertex(evt);
        }
    }
}
