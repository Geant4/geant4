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
// $Id: G4SamplingPostStepAction.cc 66241 2012-12-13 18:34:42Z gunter $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4SamplingPostStepAction.cc
//
// ----------------------------------------------------------------------

#include "G4SamplingPostStepAction.hh"
#include "G4Track.hh"
#include "G4ParticleChange.hh"
#include "G4VImportanceSplitExaminer.hh"
#include "G4Nsplit_Weight.hh"
#include "G4VTrackTerminator.hh"
#include <sstream>

G4SamplingPostStepAction::
G4SamplingPostStepAction(const G4VTrackTerminator &TrackTerminator)
  : fTrackTerminator(TrackTerminator)
{
}

G4SamplingPostStepAction::~G4SamplingPostStepAction()
{
}

void G4SamplingPostStepAction::DoIt(const G4Track& aTrack, 
                                          G4ParticleChange *aParticleChange,
                                    const G4Nsplit_Weight &nw)
{  
  // evaluate results from sampler
  if (nw.fN>1)
  {
    // split track 
    Split(aTrack, nw, aParticleChange);
  }
  else if (nw.fN==1)
  {
    // don't split, but weight may be changed ! 
    aParticleChange->ProposeWeight(nw.fW);
  }
  else if (nw.fN==0)
  {
    // kill track
    fTrackTerminator.KillTrack();
  }
  else
  {
    // wrong answer
    std::ostringstream os;
    os << "Sampler returned nw = "
       << nw
       << "\n";
    G4String msg = os.str();
    
    G4Exception("G4SamplingPostStepAction::DoIt()",
                "InvalidCondition", FatalException, msg);
  }
}

void G4SamplingPostStepAction::Split(const G4Track &aTrack,
                                     const G4Nsplit_Weight &nw,
                                           G4ParticleChange *aParticleChange)
{
  aParticleChange->ProposeWeight(nw.fW);
  aParticleChange->SetNumberOfSecondaries(nw.fN-1);
  
  for (G4int i=1;i<nw.fN;i++)
  {
    G4Track *ptrack = new G4Track(aTrack);
    
    //    ptrack->SetCreatorProcess(aTrack.GetCreatorProcess());
    ptrack->SetWeight(nw.fW);
    
    if (ptrack->GetMomentumDirection() != aTrack.GetMomentumDirection())
    {
      G4Exception("G4SamplingPostStepAction::Split()", "InvalidCondition",
                  FatalException, "Track with same momentum !");
    }
    aParticleChange->AddSecondary(ptrack);
  }
  return;
}  
