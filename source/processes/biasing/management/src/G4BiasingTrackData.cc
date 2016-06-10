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
#include "G4BiasingTrackData.hh"
#include "G4BiasingTrackDataStore.hh"

G4BiasingTrackData::G4BiasingTrackData(const G4Track* track)
  : fTrack(track),
    fBirthOperation (0),
    fBirthOperator  (0),
    fBirthProcess(0)
{
  G4BiasingTrackDataStore::GetInstance()->Register( this );
}

G4BiasingTrackData::G4BiasingTrackData(const G4Track*                   track,
				       const G4VBiasingOperation*       birthOperation,
				       const G4VBiasingOperator*        birthOperator,
				       const G4BiasingProcessInterface* birthProcess)
  : fTrack(track),
    fBirthOperation (birthOperation),
    fBirthOperator  (birthOperator),
    fBirthProcess(birthProcess)
{
  G4BiasingTrackDataStore::GetInstance()->Register( this );
}

G4BiasingTrackData::~G4BiasingTrackData()
{
  G4BiasingTrackDataStore::GetInstance()->DeRegister( this );
}


