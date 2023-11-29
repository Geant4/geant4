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
// ------------------------------------------------------------
//      GEANT 4 class implementation file
// ------------------------------------------------------------
//

#include "G4ErrorTrackLengthTarget.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"

#ifdef G4VERBOSE
#  include "G4ErrorPropagatorData.hh"  //for verbosity checking
#endif

//----------------------------------------------------------------------------
G4ErrorTrackLengthTarget::G4ErrorTrackLengthTarget(const G4double maxTrkLength)
  : G4VDiscreteProcess("G4ErrorTrackLengthTarget")
  , theMaximumTrackLength(maxTrkLength)
{
  theType = G4ErrorTarget_TrkL;

  G4ParticleTable::G4PTblDicIterator* theParticleIterator =
    G4ParticleTable::GetParticleTable()->GetIterator();

  // loop over all particles in G4ParticleTable

  theParticleIterator->reset();
  while((*theParticleIterator)())  // Loop checking, 06.08.2015, G.Cosmo
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager     = particle->GetProcessManager();
    if(!particle->IsShortLived())
    {
      // Add transportation process for all particles other than  "shortlived"
      if(pmanager == 0)
      {
        // Error !! no process manager
        G4String particleName = particle->GetParticleName();
        G4Exception("G4ErrorTrackLengthTarget::G4ErrorTrackLengthTarget",
                    "No process manager", RunMustBeAborted, particleName);
      }
      else
      {
        G4ProcessVector* procvec = pmanager->GetProcessList();
        std::size_t isiz = (G4int)procvec->size();

        for(G4int ii = 0; ii < (G4int)isiz; ++ii)
        {
          if(((*procvec)[ii])->GetProcessName() == "G4ErrorTrackLengthTarget")
          {
            pmanager->RemoveProcess((*procvec)[ii]);
          }
        }
        pmanager->AddDiscreteProcess(this, 4);
        isiz = procvec->size();
      }
    }
    else
    {
      // shortlived particle case
    }
  }
}

//-----------------------------------------------------------------------
G4double G4ErrorTrackLengthTarget::PostStepGetPhysicalInteractionLength(
  const G4Track& track, G4double, G4ForceCondition* condition)
{
  *condition = NotForced;
  return GetMeanFreePath(track, 0., condition);
}

//-----------------------------------------------------------------------
G4double G4ErrorTrackLengthTarget::GetMeanFreePath(const class G4Track& track,
                                                   G4double,
                                                   enum G4ForceCondition*)
{
#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 3)
  {
    G4cout << " G4ErrorTrackLengthTarget::GetMeanFreePath "
           << theMaximumTrackLength - track.GetTrackLength() << G4endl;
  }
#endif

  return theMaximumTrackLength - track.GetTrackLength();
}

G4VParticleChange* G4ErrorTrackLengthTarget::PostStepDoIt(const G4Track& aTrack,
                                                          const G4Step&)
{
  theParticleChange.Initialize(aTrack);
  return &theParticleChange;
}

//-----------------------------------------------------------------------
void G4ErrorTrackLengthTarget::Dump(const G4String& msg) const
{
  G4cout << msg << "G4ErrorTrackLengthTarget: max track length = "
         << theMaximumTrackLength << G4endl;
}
