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
// $Id: G4ParallelScoreSampler.cc,v 1.4 2002-09-02 13:27:26 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelScoreSampler.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelScoreSampler.hh"
#include "G4PScoreProcess.hh"
#include "G4ProcessPlacer.hh"
#include "G4ParallelWorld.hh"
#include "G4ParallelTransport.hh"


G4ParallelScoreSampler::
G4ParallelScoreSampler(G4VPhysicalVolume &worldvolume,
		       const G4String &particlename,
		       G4VPScorer &scorer) : 
  fParticleName(particlename),
  fParallelWorld(new G4ParallelWorld(worldvolume)),
  fPScorer(scorer),
  fPScorerProcess(0),
  fParallelTransport(0)
{}

G4ParallelScoreSampler::~G4ParallelScoreSampler()
{
  if (fPScorerProcess) {
    G4ProcessPlacer placer(fParticleName);
    placer.RemoveProcess(fPScorerProcess);
    delete  fPScorerProcess;
  }
  delete fParallelWorld;
}


G4PScoreProcess *G4ParallelScoreSampler::CreateParallelScoreProcess()
{
  if (!fPScorerProcess) {
    fPScorerProcess = 
      new G4PScoreProcess(fParallelWorld->GetParallelStepper(), 
			  fPScorer);
  }
  return fPScorerProcess;
}

G4ParallelTransport *G4ParallelScoreSampler::CreateParallelTransport()
{
  if (!fParallelTransport) {
    fParallelTransport = new G4ParallelTransport(fParallelWorld->GetGeoDriver(), 
						 fParallelWorld->
						 GetParallelStepper());
  }
  return fParallelTransport;
}

void G4ParallelScoreSampler::Initialize()
{
  G4ProcessPlacer placer(fParticleName);
  placer.AddProcessAsSecondDoIt(CreateParallelScoreProcess());
  placer.AddProcessAsSecondDoIt(CreateParallelTransport());
}




