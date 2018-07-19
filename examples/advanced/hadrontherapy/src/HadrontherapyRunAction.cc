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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#include "HadrontherapyRunAction.hh"
#include "HadrontherapyEventAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4Timer.hh"
#include "HadrontherapyRunAction.hh"

#include "HadrontherapyMatrix.hh"

HadrontherapyRunAction::HadrontherapyRunAction()
{

    
}

HadrontherapyRunAction::~HadrontherapyRunAction()
{
    
}



void HadrontherapyRunAction::BeginOfRunAction(const G4Run* aRun)
{
    G4RunManager::GetRunManager()-> SetRandomNumberStore(true);
    G4cout << "Run " << aRun -> GetRunID() << " starts ..." << G4endl;
    
    electromagnetic = 0;
    hadronic = 0;
    
    
    
}

void HadrontherapyRunAction::EndOfRunAction(const G4Run*)
{
    //G4cout << " Summary of Run " << aRun -> GetRunID() <<" :"<< G4endl;
    //G4cout << "Number of electromagnetic processes of primary particles in the phantom:"
    // 	   << electromagnetic << G4endl;
    //G4cout << "Number of hadronic processes of primary particles in the phantom:"
    //	   << hadronic << G4endl;
}
void HadrontherapyRunAction::AddEMProcess()
{
    electromagnetic += 1;
}
void HadrontherapyRunAction::AddHadronicProcess()
{
    hadronic += 1;
}

