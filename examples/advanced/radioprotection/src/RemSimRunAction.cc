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
//    *******************************
//    *                             *
//    *    RemSimRunAction.cc       *
//    *                             *
//    *******************************
//
// Code developed by: S.Guatelli, guatelli@ge.infn.it
// $Id: RemSimRunAction.cc,v 1.14 2006-06-29 16:24:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "RemSimRunAction.hh"
#include "RemSimDetectorConstruction.hh"
#include "RemSimPrimaryGeneratorAction.hh"
#ifdef G4ANALYSIS_USE
#include "RemSimAnalysisManager.hh"
#endif
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "RemSimRunAction.hh"

RemSimRunAction::RemSimRunAction()
{
}

RemSimRunAction::~RemSimRunAction()
{
 
 }

void RemSimRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun -> GetRunID() << " start." << G4endl;
#ifdef G4ANALYSIS_USE
  G4int runNb = aRun -> GetRunID();
  if (runNb == 0) 
    {
     RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
     analysis -> book();
    }
#endif
}

void RemSimRunAction::EndOfRunAction(const G4Run* aRun)
{  
 G4double numberEvents = aRun -> GetNumberOfEvent();
 G4cout<< "Number of events:" << numberEvents << G4endl;
}

