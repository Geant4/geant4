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
//
// $Id: G4ImportancePostStepDoIt.cc,v 1.8 2002-11-04 10:47:56 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ImportancePostStepDoIt.cc
//
// ----------------------------------------------------------------------

#include "G4ImportancePostStepDoIt.hh"
#include "G4Track.hh"
#include "G4ParticleChange.hh"
#include "G4VImportanceSplitExaminer.hh"
#include "G4Nsplit_Weight.hh"
#include "G4VTrackTerminator.hh"
#include "g4std/strstream"


G4ImportancePostStepDoIt::
G4ImportancePostStepDoIt(const G4VTrackTerminator &TrackTerminator)
  :
  fTrackTerminator(TrackTerminator)
{}
G4ImportancePostStepDoIt::~G4ImportancePostStepDoIt(){}

void G4ImportancePostStepDoIt::DoIt(const G4Track& aTrack, 
				    G4ParticleChange *aParticleChange,
				    const G4Nsplit_Weight &nw)
{  
  // evaluate results from sampler
  if (nw.fN>1) {
    // split track 
    Split(aTrack, nw, aParticleChange);
  }
  else if (nw.fN==1) {
    // don't split, but weight may be changed ! 
    aParticleChange->SetWeightChange(nw.fW);
  }
  else if (nw.fN==0) {
    // kill track
    fTrackTerminator.KillTrack();
  }
  else {
    // wrong answer
    char st[200];
    G4std::ostrstream os(st,200);
    os << "G4ImportancePostStepDoIt::DoIt: sampler returned nw = "
       << nw
       << "\n"
       << '\0';
    G4String m(st);
    
    G4Exception(m);
  }
}

void G4ImportancePostStepDoIt::Split(const G4Track &aTrack,
				     const G4Nsplit_Weight &nw,
				     G4ParticleChange *aParticleChange)
{
  aParticleChange->SetWeightChange(nw.fW);
  aParticleChange->SetNumberOfSecondaries(nw.fN-1);
  
  for (G4int i=1;i<nw.fN;i++) {
    G4Track *ptrack = new G4Track(aTrack);
    
    //    ptrack->SetCreatorProcess(aTrack.GetCreatorProcess());
    ptrack->SetWeight(nw.fW);
    
    if (ptrack->GetMomentumDirection() != aTrack.GetMomentumDirection()) {
      G4Exception("ERROR - G4ImportancePostStepDoIt::Split: (ptrack->GetMomentumDirection() != aTrack.GetMomentumDirection()");
    }
    
    aParticleChange->AddSecondary(ptrack);
  }
  return;
}  
