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
// $Id: BrachyRunAction.cc,v 1.18 2006-06-29 15:48:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "BrachyRunAction.hh"

#ifdef G4ANALYSIS_USE
#include "BrachyAnalysisManager.hh"
#endif

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "BrachyRunAction.hh"

BrachyRunAction::BrachyRunAction()
{
}

BrachyRunAction::~BrachyRunAction()
{ 
}
void BrachyRunAction::BeginOfRunAction(const G4Run* aRun)
{ 
 G4cout << "### Run " << aRun -> GetRunID() << " start." << G4endl;

#ifdef G4ANALYSIS_USE
 G4int runNb = aRun -> GetRunID();
 if (runNb == 0) 
    {  
     BrachyAnalysisManager* analysis = BrachyAnalysisManager::getInstance();
     analysis->book();
    }
 else { G4cout << "The results of Run:"<< runNb << " are summed to the" << 
        " results of the previous Run in brachytherapy.hbk" << G4endl;} 
#endif  
}

void BrachyRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout << "number of events = " << aRun->GetNumberOfEvent() << G4endl;
}




