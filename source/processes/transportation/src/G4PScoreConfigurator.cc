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
// $Id: G4PScoreConfigurator.cc,v 1.9 2006/11/14 09:11:18 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-02 $
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
  
