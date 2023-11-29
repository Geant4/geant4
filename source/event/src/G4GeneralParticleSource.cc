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
// G4GeneralParticleSource class implementation
//
// Author: Fan Lei, QinetiQ ltd - 05/02/2004
// Customer: ESA/ESTEC
// Version: 2.0
// --------------------------------------------------------------------

#include "G4Event.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"
#include "G4SingleParticleSource.hh"
#include "G4UnitsTable.hh"

#include "G4GeneralParticleSourceData.hh"

#include "G4Threading.hh"
#include "G4AutoLock.hh"

namespace
{
  G4Mutex messangerInit = G4MUTEX_INITIALIZER;
}

G4GeneralParticleSource::G4GeneralParticleSource()
{
  GPSData = G4GeneralParticleSourceData::Instance();
  // currentSource = GPSData->GetCurrentSource();
  // currentSourceIdx = G4int(GPSData->GetSourceVectorSize() - 1);

  // Messenger is special, only a worker should instantiate it.
  // Singleton pattern
  //
  theMessenger = G4GeneralParticleSourceMessenger::GetInstance(this);

  // Some initialization should be done only once
  //
  G4AutoLock l(&messangerInit);
  static G4bool onlyOnce = false;
  if ( !onlyOnce )
  {
    theMessenger->SetParticleGun(GPSData->GetCurrentSource());
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
  GPSData->Lock();

  GPSData->AddASource(aV);
  theMessenger->SetParticleGun(GPSData->GetCurrentSource());

  // TODO: But do we really normalize here after each source?
  IntensityNormalization();

  GPSData->Unlock();
}

void G4GeneralParticleSource::IntensityNormalization()
{
  GPSData->IntensityNormalise();
  normalised=GPSData->Normalised();
}

void G4GeneralParticleSource::ListSource()
{
  G4cout << "The number of particle sources is: "
         << GPSData->GetIntensityVectorSize() << G4endl;
  G4cout << " Multiple Vertex sources: " << GPSData->GetMultipleVertex();
  G4cout << " Flat Sampling flag: " << GPSData->GetFlatSampling() << G4endl;
  const G4int currentIdx = GPSData->GetCurrentSourceIdx();
  for(G4int i=0; i<GPSData->GetIntensityVectorSize(); ++i)
  {
    G4cout << "\tsource " << i << " with intensity: "
           << GPSData->GetIntensity(i) << G4endl;
    const G4SingleParticleSource* thisSrc = GPSData->GetCurrentSource(i);
    G4cout << " \t\tNum Particles: "<<thisSrc->GetNumberOfParticles()
           << "; Particle type: "
           << thisSrc->GetParticleDefinition()->GetParticleName() << G4endl;
    G4cout << " \t\tEnergy: "
           << G4BestUnit(thisSrc->GetParticleEnergy(),"Energy") << G4endl;
    G4cout << " \t\tDirection: "
           << thisSrc->GetAngDist()->GetDirection() << "; Position: ";
    G4cout << G4BestUnit(thisSrc->GetPosDist()->GetCentreCoords(),"Length")
           << G4endl;
    G4cout << " \t\tAngular Distribution: "
           << thisSrc->GetAngDist()->GetDistType() << G4endl;
    G4cout << " \t\tEnergy Distribution: "
           << thisSrc->GetEneDist()->GetEnergyDisType() << G4endl;
    G4cout << " \t\tPosition Distribution Type: "
           << thisSrc->GetPosDist()->GetPosDisType();
    G4cout << "; Position Shape: "
           << thisSrc->GetPosDist()->GetPosDisShape() << G4endl;
  }

  // Set back previous source
  GPSData->GetCurrentSource(currentIdx);
}

void G4GeneralParticleSource::SetCurrentSourceto(G4int aV)
{
  G4int id = aV;
  if ( id < GPSData->GetIntensityVectorSize() )
  {
    // currentSourceIdx = aV;
    // currentSource = GPSData->GetCurrentSource(id);
    theMessenger->SetParticleGun(GPSData->GetCurrentSource(id));
  }
  else
  {
    G4ExceptionDescription msg;
    msg << "Trying to set source to index " << aV << " but only "
        << GPSData->GetIntensityVectorSize() << " sources are defined.";
    G4Exception("G4GeneralParticleSoruce::SetCurrentSourceto", "G4GPS004",
                FatalException, msg);
  }
}

void G4GeneralParticleSource::SetCurrentSourceIntensity(G4double aV)
{
  GPSData->Lock();
  GPSData->SetCurrentSourceIntensity(aV);
  GPSData->Unlock();
  normalised = GPSData->Normalised();
}

void G4GeneralParticleSource::ClearAll()
{
  GPSData->ClearSources();
  normalised=GPSData->Normalised();
}

void G4GeneralParticleSource::DeleteaSource(G4int aV)
{
  G4int id = aV;
  if ( id <= GPSData->GetIntensityVectorSize() )
  {
    GPSData->DeleteASource(aV);
    normalised=GPSData->Normalised();
  }
  else
  {
    G4cout << " source index is invalid " << G4endl;
    G4cout << "    it shall be <= "
           << GPSData->GetIntensityVectorSize() << G4endl;
  }
}

void G4GeneralParticleSource::GeneratePrimaryVertex(G4Event* evt)
{
  if (!GPSData->GetMultipleVertex())
  {
    G4SingleParticleSource* currentSource = GPSData->GetCurrentSource();
    if (GPSData->GetIntensityVectorSize() > 1)
    {
      // Try to minimize locks
      if (! normalised )
      {
        // According to local variable, normalization is needed
        // Check with underlying shared resource, another
        // thread could have already normalized this
        GPSData->Lock();
        G4bool norm = GPSData->Normalised();
        if (!norm)
        {
          IntensityNormalization();
        }
        // This takes care of the case in which the local variable
        // is False and the underlying resource is true.
        normalised = GPSData->Normalised();
        GPSData->Unlock();
      }
      G4double rndm = G4UniformRand();
      G4int i = 0 ;
      if (! GPSData->GetFlatSampling() )
      {
        while ( rndm > GPSData->GetSourceProbability(i) ) ++i;
          currentSource = GPSData->GetCurrentSource(i);
      }
      else
      {
        i = G4int (GPSData->GetIntensityVectorSize()*rndm);
        currentSource = GPSData->GetCurrentSource(i);
      }
    }
    currentSource->GeneratePrimaryVertex(evt);
  } 
  else
  {
    for (G4int i = 0; i <  GPSData->GetIntensityVectorSize(); ++i)
    {
      GPSData->GetCurrentSource(i)->GeneratePrimaryVertex(evt);
    }
  }
}
