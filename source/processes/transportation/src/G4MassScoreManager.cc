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
// $Id: G4MassScoreManager.cc,v 1.2 2002-04-09 17:40:15 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4MassScoreManager.cc
//
// ----------------------------------------------------------------------

#include "G4MassScoreManager.hh"
#include "G4MScoreProcess.hh"
#include "G4ProcessPlacer.hh"


G4MassScoreManager::G4MassScoreManager(G4VPScorer &ascorer,
                                       const G4String &particlename)
 : fScorer(ascorer),
   fParticleName(particlename),
   fMScoreProcess(0)
{}

G4MassScoreManager::~G4MassScoreManager()
{
  if (fMScoreProcess) delete fMScoreProcess;
}

G4MScoreProcess *G4MassScoreManager::CreateMassScoreProcess()
{
  if (!fMScoreProcess) {
    fMScoreProcess = new G4MScoreProcess(fScorer);
  }
  return fMScoreProcess;
}

void G4MassScoreManager::Initialize()
{
  G4ProcessPlacer placer(fParticleName);
  placer.AddProcessAsSecondDoIt(CreateMassScoreProcess());
}
