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
// $Id: G4ScoreConfigurator.cc,v 1.1 2006/11/20 10:02:17 ahoward Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// Class G4ScoreConfigurator
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#include "G4ScoreConfigurator.hh"
#include "G4VScorer.hh"
#include "G4VTrackTerminator.hh"
#include "G4ScoreProcess.hh"

G4ScoreConfigurator::
G4ScoreConfigurator(G4VPhysicalVolume* worldvolume, const G4String &particlename,
                           G4VScorer &scorer, G4bool para) 
  : fWorld(worldvolume),
    fPlacer(particlename),
    fScorer(scorer),
    fScoreProcess(0)
{
  paraflag = para;
}

G4ScoreConfigurator::~G4ScoreConfigurator()
{
  if (fScoreProcess)
  {
    fPlacer.RemoveProcess(fScoreProcess);
    delete fScoreProcess;
  }
}
  
const G4VTrackTerminator *
G4ScoreConfigurator::GetTrackTerminator() const
{
  return fScoreProcess;
}

void G4ScoreConfigurator::Configure(G4VSamplerConfigurator *)
{
  G4cout << " entering scoring configure " << G4endl;
  fScoreProcess = new G4ScoreProcess(fScorer,"ScoringProcess",paraflag);
  G4cout << " creating score process " << G4endl;
  if (!fScoreProcess)
  {
    G4Exception("G4ScoreConfigurator::Configure()",
                "FatalError", FatalException,
                "Failed allocation of G4ScoreProcess !");
  }

  G4cout << " setting parallel world " << G4endl;
  if(paraflag) fScoreProcess->SetParallelWorld(fWorld);
  G4cout << " set " << G4endl;

  fPlacer.AddProcessAsSecondDoIt(fScoreProcess);
}
