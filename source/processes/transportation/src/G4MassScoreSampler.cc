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
// $Id: G4MassScoreSampler.cc,v 1.1 2002-05-31 10:16:02 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4MassScoreSampler.cc
//
// ----------------------------------------------------------------------

#include "G4MassScoreSampler.hh"
#include "G4MScoreProcess.hh"
#include "G4ProcessPlacer.hh"


G4MassScoreSampler::G4MassScoreSampler(G4VPScorer &ascorer,
                                       const G4String &particlename)
 : fScorer(ascorer),
   fParticleName(particlename),
   fMScoreProcess(0)
{}

G4MassScoreSampler::~G4MassScoreSampler()
{
  if (fMScoreProcess) {
    G4ProcessPlacer placer(fParticleName);
    placer.RemoveProcess(fMScoreProcess);
    delete fMScoreProcess;
  }
}

G4MScoreProcess *G4MassScoreSampler::CreateMassScoreProcess()
{
  if (!fMScoreProcess) {
    fMScoreProcess = new G4MScoreProcess(fScorer);
  }
  return fMScoreProcess;
}

void G4MassScoreSampler::Initialize()
{
  G4ProcessPlacer placer(fParticleName);
  placer.AddProcessAsSecondDoIt(CreateMassScoreProcess());
}
