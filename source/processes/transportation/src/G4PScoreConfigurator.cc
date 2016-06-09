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
// $Id: G4PScoreConfigurator.cc,v 1.7 2003/11/26 14:51:49 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// ----------------------------------------------------------------------
// Class G4PScoreConfigurator
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#include "G4PScoreConfigurator.hh"

#include "G4VTrackTerminator.hh"
#include "G4PScoreProcess.hh"

G4PScoreConfigurator::
G4PScoreConfigurator(const G4String &particlename,
                           G4VParallelStepper &pstepper,
                           G4VScorer &scorer) 
  : fPlacer(particlename),
    fPStepper(pstepper),
    fScorer(scorer),
    fPScoreProcess(0)
{
}

G4PScoreConfigurator::~G4PScoreConfigurator()
{
  if (fPScoreProcess)
  {
    fPlacer.RemoveProcess(fPScoreProcess);
    delete fPScoreProcess;
  }
}
 
void G4PScoreConfigurator::Configure(G4VSamplerConfigurator *)
{
  fPScoreProcess = new G4PScoreProcess(fPStepper, fScorer);
  if (!fPScoreProcess)
  {
    G4Exception("G4PScoreConfigurator::Configure()", "FatalError",
                FatalException, "Failed to allocate G4PScoreProcess !");
  }
  fPlacer.AddProcessAsSecondDoIt(fPScoreProcess);
}

const G4VTrackTerminator *
G4PScoreConfigurator::GetTrackTerminator() const
{
  return fPScoreProcess;
}
  
