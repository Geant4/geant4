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
// $Id: BrachyRunAction.cc,v 1.2 2004/05/25 08:36:18 guatelli Exp $
// GEANT4 tag $Name: geant4-07-01 $
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




