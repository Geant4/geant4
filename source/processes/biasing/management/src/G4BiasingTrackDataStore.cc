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
#include "G4BiasingTrackDataStore.hh"
#include "G4BiasingTrackData.hh"


//G4BiasingTrackDataStore* G4BiasingTrackDataStore::fInstance = 0;

G4BiasingTrackDataStore* G4BiasingTrackDataStore::GetInstance()
{
    static G4ThreadLocalSingleton<G4BiasingTrackDataStore> instance;
    return instance.Instance();
}

void G4BiasingTrackDataStore::Register(G4BiasingTrackData* data)
{
  fTrackDataStore[data->GetTrack()] = data;
}

void G4BiasingTrackDataStore::DeRegister(G4BiasingTrackData* data)
{
  fTrackDataStore[data->GetTrack()] = 0;
}

G4BiasingTrackDataStore::G4BiasingTrackDataStore()
{}

G4BiasingTrackDataStore::~G4BiasingTrackDataStore()
{
  for ( std::map < const G4Track*, G4BiasingTrackData* >::iterator it = fTrackDataStore.begin() ;
	it != fTrackDataStore.end() ; it++ )
    {
      G4BiasingTrackData* data = (*it).second;
      if ( data != 0 ) delete data;
    }
}

