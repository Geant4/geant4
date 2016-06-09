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
// $Id: G4MScoreConfigurator.cc,v 1.9 2006/11/13 16:17:16 japost Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// Class G4MScoreConfigurator
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#include "G4MScoreConfigurator.hh"
#include "G4VScorer.hh"
#include "G4VTrackTerminator.hh"
#include "G4MScoreProcess.hh"

G4MScoreConfigurator::
G4MScoreConfigurator(const G4String &particlename,
                           G4VScorer &scorer) 
  : fPlacer(particlename),
    fScorer(scorer),
    fMScoreProcess(0)
{
}

G4MScoreConfigurator::~G4MScoreConfigurator()
{
  if (fMScoreProcess)
  {
    fPlacer.RemoveProcess(fMScoreProcess);
    delete fMScoreProcess;
  }
}
  
const G4VTrackTerminator *
G4MScoreConfigurator::GetTrackTerminator() const
{
  return fMScoreProcess;
}

void G4MScoreConfigurator::Configure(G4VSamplerConfigurator *)
{
  fMScoreProcess = new G4MScoreProcess(fScorer);
  if (!fMScoreProcess)
  {
    G4Exception("G4MScoreConfigurator::Configure()",
                "FatalError", FatalException,
                "Failed allocation of G4MScoreProcess !");
  }
  fPlacer.AddProcessAsSecondDoIt(fMScoreProcess);
}
