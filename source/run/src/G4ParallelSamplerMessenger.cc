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
// $Id: G4ParallelSamplerMessenger.cc,v 1.3 2002-07-29 16:03:15 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelSamplerMessenger.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelSamplerMessenger.hh"
#include "g4std/fstream"

#include "G4ImportanceGeometryConstructor.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParallelImportanceScoreSampler.hh"
#include "G4ParallelScoreSampler.hh"
#include "G4VPrintableScorer.hh"
#include "G4VPrintableScorerFactory.hh"

G4ParallelSamplerMessenger::
G4ParallelSamplerMessenger(G4ImportanceGeometryConstructor *igc,
			   const G4VPrintableScorerFactory &psf) : 
  fImpGeoConst(igc),
  fPrintableScorerFactory(psf) 
{
  fImpCmd = new G4UIcmdWithAString("/imp/importancesample", this);
  fScoreCmd = new G4UIcmdWithAString("/imp/score", this);
  fPrintCmd = new G4UIcmdWithAString("/imp/writefile", this);
  fImpSampler = 0;
  fImpPartilce = "none";
}

void G4ParallelSamplerMessenger::
SetNewValue(G4UIcommand * command, G4String newValue){
  if (command==fImpCmd) {
    ImpCmdAction(newValue);
  }
  if (command==fScoreCmd) {
    ScoreCmdAction(newValue);
  }
  if (command==fPrintCmd) {
    PrintCmdAction(newValue);
  }
};
  
void G4ParallelSamplerMessenger::
ImpCmdAction(const G4String &particlename) {
  if (fImpSampler!=0) {
    Error(G4String("SetNewValue: an importance sampler already exists\n")
	  + 
	  G4String("only one importance sampler can")
	  + G4String("be used with this messenger"));
  }
  CheckNewParticle(particlename);
  fImpPartilce = particlename;
  fImpScorer = fPrintableScorerFactory.
    CreatePrintableScorer(particlename);
  fImpSampler = new 
    G4ParallelImportanceScoreSampler(*(fImpGeoConst->GetIStore()),
				     *fImpScorer->GetPointerToScorer(),
				     particlename);
  fMapNameSampler[particlename] = fImpSampler;
  fImpSampler->Initialize();
  
}
void G4ParallelSamplerMessenger::
ScoreCmdAction(const G4String &particlename){
  CheckNewParticle(particlename);
  G4VPrintableScorer *s = 
    fPrintableScorerFactory.
    CreatePrintableScorer(particlename);
  fMapNameScorer[particlename] = s;
  fMapNameSampler[particlename] = new 
    G4ParallelScoreSampler(*(fImpGeoConst->GetWorldVolume()),
			   particlename,
			   *s->GetPointerToScorer());
  fMapNameSampler[particlename]->Initialize();
}
void G4ParallelSamplerMessenger::
PrintCmdAction(const G4String &filename){
  G4std::ofstream of(filename);
  if (fImpSampler) {
    of << "  Results for importanct sampled particle: "
       <<  fImpPartilce << G4endl;
    of << "  --------------------------------------------------"<<G4endl; 
    fImpScorer->PrintTable(&of, fImpGeoConst->GetIStore());
    of << "----------------------------------------------------" 
       << "----------------------------------------------------"
       << "---------------------"
       << G4endl;
    of << G4endl;
  }
  for (G4MapNameScorer::iterator it = fMapNameScorer.begin();
       it != fMapNameScorer.end(); it++) {
    of << "  Results for scored particle: " 
       << it->first << G4endl;
    of << "  --------------------------------------------------"<<G4endl; 
    it->second->PrintTable(&of);
    of << "----------------------------------------------------" 
       << "----------------------------------------------------"
       << "---------------------"
       << G4endl;
    of << G4endl;
  }
}

void G4ParallelSamplerMessenger::
CheckNewParticle(const G4String &particlename){
  G4MapNameSampler::iterator it = fMapNameSampler.find(particlename);
  if (it!=fMapNameSampler.end()) {
    Error("CheckNewParticle: particle:" + particlename + 
	  ", already booked for sampling");
  }
};
