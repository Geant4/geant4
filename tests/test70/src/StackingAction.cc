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
#include "G4TrackStatus.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "StackingAction.hh"
#include "G4EventManager.hh"
#include "G4ITManager.hh"
#include "G4ITStepManager.hh"

StackingAction::StackingAction(G4bool flag)
{
    fLaunchITStepMan = flag;
}


StackingAction::~StackingAction()
{;}

void StackingAction::NewStage()
{
    if(stackManager->GetNTotalTrack() == 0 && fLaunchITStepMan )
    {
        G4cout << "You are launching the chemistry ..." << G4endl ;
        G4ITStepManager::Instance() -> Process();
        G4ITStepManager::Instance() -> ClearList();
    }
}
