//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

#include "ExN07StackingAction.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

G4int ExN07StackingAction::nGamma[6] = {0,0,0,0,0,0};
G4int ExN07StackingAction::nElectron[6] = {0,0,0,0,0,0};
G4int ExN07StackingAction::nPositron[6] = {0,0,0,0,0,0};
G4double ExN07StackingAction::eMinGamma[6] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
G4double ExN07StackingAction::eMinElectron[6] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};
G4double ExN07StackingAction::eMinPositron[6] = {DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};

ExN07StackingAction::ExN07StackingAction()
{;}

ExN07StackingAction::~ExN07StackingAction()
{;}

G4ClassificationOfNewTrack 
ExN07StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4TouchableHistory* theTouchable
   = (G4TouchableHistory*)(aTrack->GetTouchable());
  if(theTouchable&&theTouchable->GetHistoryDepth()>1)
  {
    G4VPhysicalVolume* calPhys = theTouchable->GetVolume(1);
    if(calPhys)
    {
      G4int i = 0;
      if(calPhys->GetName()=="Cal-B") { i = 2; }
      else if(calPhys->GetName()=="Cal-C") { i = 4; }
      G4int layerNumber = theTouchable->GetReplicaNumber();
      i += layerNumber%2;
      G4ParticleDefinition * particleType = aTrack->GetDefinition();
      if(particleType==G4Gamma::GammaDefinition())
      {
        nGamma[i]++;
        if(eMinGamma[i]>aTrack->GetKineticEnergy())
        { eMinGamma[i]=aTrack->GetKineticEnergy(); }
      }
      else if(particleType==G4Electron::ElectronDefinition())
      {
        nElectron[i]++;
        if(eMinElectron[i]>aTrack->GetKineticEnergy())
        { eMinElectron[i]=aTrack->GetKineticEnergy(); }
      }
      else if(particleType==G4Positron::PositronDefinition())
      {
        nPositron[i]++;
        if(eMinPositron[i]>aTrack->GetKineticEnergy())
        { eMinPositron[i]=aTrack->GetKineticEnergy(); }
      }
    }
  }
  return fUrgent;
}

void ExN07StackingAction::PrepareNewEvent()
{ 
  for(int i=0;i<6;i++)
  {
    nGamma[i] = 0; 
    nElectron[i] = 0;
    nPositron[i] = 0;
    eMinGamma[i] = DBL_MAX;
    eMinElectron[i] = DBL_MAX;
    eMinPositron[i] = DBL_MAX;
  }
}


