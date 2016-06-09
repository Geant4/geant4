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
// --------------------------------------------------------------
//                 GEANT 4 - Brachytherapy example
// --------------------------------------------------------------
//
// Code developed by:
//  S.Guatelli
//
//
//    *******************************
//    *                             *
//    *    BrachyRunAction.cc       *
//    *                             *
//    *******************************
//
// $Id: BrachyRunAction.cc,v 1.3 2006/06/29 17:33:38 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//

#include "BrachyRunAction.hh"
#include "BrachyEventAction.hh"

#ifdef G4ANALYSIS_USE
#include "BrachyAnalysisManager.hh"
#endif

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "BrachyDetectorConstruction.hh"
#include "BrachyRunMessenger.hh"
#include "G4SDManager.hh"
#include "G4Timer.hh"
#include "BrachyFactoryIr.hh"
#include "BrachyFactoryI.hh"
#include "BrachyFactory.hh"
#include "BrachyRunAction.hh"

BrachyRunAction::BrachyRunAction()
{
  runMessenger = new BrachyRunMessenger(this);
}

BrachyRunAction::~BrachyRunAction()
{ 
  delete runMessenger;
}
void BrachyRunAction::BeginOfRunAction(const G4Run*)
{ 
#ifdef G4ANALYSIS_USE
  BrachyAnalysisManager* analysis = BrachyAnalysisManager::getInstance();
  analysis -> book();
#endif  
  G4RunManager* runManager = G4RunManager::GetRunManager();

  if(runManager)
    { switch(sourceChoice)
      {
        case 1:
	  factory = new BrachyFactoryI;
	  break;
        default:   
	  factory = new BrachyFactoryIr; 
      }      
    G4VUserPrimaryGeneratorAction* sourcePrimaryParicle = 
                                     factory -> CreatePrimaryGeneratorAction();
      
    if(sourcePrimaryParicle) runManager -> SetUserAction(sourcePrimaryParicle);     
    }
}

void BrachyRunAction::SelectEnergy(G4int choice)
{
  sourceChoice = choice;
  if (sourceChoice == 1) factory = new BrachyFactoryI;
  else factory = new BrachyFactoryIr; 
}

void BrachyRunAction::EndOfRunAction(const G4Run* aRun)
{
#ifdef G4ANALYSIS_USE
  BrachyAnalysisManager* analysis = BrachyAnalysisManager::getInstance();
#endif
  G4cout << "number of event = " << aRun -> GetNumberOfEvent() << G4endl;
  
#ifdef G4ANALYSIS_USE      
  analysis -> finish();
#endif

  delete factory;    
}




