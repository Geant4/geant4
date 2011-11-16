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
// $Id: Tst33AppStarterMessenger.cc,v 1.12 2008-04-21 09:00:03 ahoward Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33AppStarterMessenger.cc
//
// ----------------------------------------------------------------------

#include "G4Types.hh"
#include <sstream>

#include "Tst33AppStarterMessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "Tst33AppStarter.hh"

Tst33AppStarterMessenger::
Tst33AppStarterMessenger(Tst33AppStarter &appstarter)
  : G4UImessenger(), fAppStarter(appstarter)
{
  fUseCoupledCmd = new G4UIcommand("/Tst33/UseCoupledTransportation", this);
  fMassGeoCmd = new G4UIcommand("/Tst33/MassGeometry", this);
  fParallelGeoCmd = new G4UIcommand("/Tst33/ParallelGeometry",this);
  fScoringCmd = new G4UIcommand("/Tst33/Scoring",this);
  fImpCmd = new G4UIcommand("/Tst33/ImportanceSampling",this);
  fWWCmd = new G4UIcmdWithAString("/Tst33/WeightWindow",this);
  fWWRCmd = new G4UIcmdWithAnInteger("/Tst33/WeightRoulette",this);
  fClearSmaplingCmd = new G4UIcommand("/Tst33/ClearSampling",this);
  fConfigureSamplingCmd = new G4UIcommand("/Tst33/ConfigureSampling",this);
  fVisAppComand = new G4UIcommand("/Tst33/Visualization",this);
  fTimedAppComand = new G4UIcmdWithAnInteger("/Tst33/Timed",this);
  fPostRunCmd = new G4UIcommand("/Tst33/PostRun",this);
  fRunCmd = new G4UIcmdWithAnInteger("/Tst33/Run", this);
  fWeightChangerCmd = new G4UIcommand("/Tst33/AddWeightChanger",this);
}

Tst33AppStarterMessenger::~Tst33AppStarterMessenger(){

  delete fUseCoupledCmd;
  delete fMassGeoCmd;
  delete fParallelGeoCmd;
  delete fScoringCmd;
  delete fImpCmd;
  delete fWWRCmd;
  delete fWWCmd;
  delete fClearSmaplingCmd;
  delete fConfigureSamplingCmd;
  delete fVisAppComand;
  delete fTimedAppComand;
  delete fPostRunCmd;
  delete fRunCmd;
  delete fWeightChangerCmd;
}

void Tst33AppStarterMessenger::SetNewValue(G4UIcommand* pCmd,
					 G4String szValue) {
  if (pCmd==fWeightChangerCmd) {
    fAppStarter.AddWeightChanger();
  }
  if (pCmd==fUseCoupledCmd) {
    fAppStarter.ForcingCoupled(true);
  }
  if (pCmd==fMassGeoCmd) {
    fAppStarter.CreateMassGeometry();
  }
  if (pCmd==fParallelGeoCmd) {
    fAppStarter.CreateParallelGeometry();
  }
  if (pCmd==fScoringCmd) {
    fAppStarter.CreateScorer();
  }
  if (pCmd==fImpCmd) {
    fAppStarter.CreateIStore();
  }
  if (pCmd==fWWCmd) {
    G4int ipoa;
    G4bool zeroWindow;
    std::istringstream is(szValue);
    is >> ipoa >>zeroWindow;

    G4PlaceOfAction poa = onBoundary;
    if (ipoa == 1) {
      poa = onBoundary;
    }
    else if (ipoa == 2) {
      poa = onCollision;
    }
    else if (ipoa == 3) {
      poa = onBoundaryAndCollision;
    }
    else {
      G4ExceptionDescription desc;
      desc << "Tst33/WeightWindow:  first argument has to be 1:"
           << " onBoundary, 2: onCollisions, 3: onBoundaryAndCollisions" << G4endl;
      desc << "first arg is : " << ipoa << ", second arg is: " << zeroWindow << G4endl;
      G4Exception("Tst33AppStarterMessenger::SetNewValue()", "TST33-02", 
               FatalException, desc);
    }
    
    fAppStarter.CreateWeightWindowStore(poa, zeroWindow);
  }
  if (pCmd==fWWRCmd) {
    G4int mode = fWWRCmd->GetNewIntValue(szValue);
    fAppStarter.CreateWeightRoulette(mode);
  }
  if (pCmd==fClearSmaplingCmd) {
    fAppStarter.ClearSampling();
  }
  if (pCmd==fConfigureSamplingCmd) {
    fAppStarter.ConfigureSampling();
  }
  if (pCmd==fVisAppComand) {
    fAppStarter.CreateVisApplication();
  }
  if (pCmd==fTimedAppComand) {
    G4int time = fTimedAppComand->GetNewIntValue(szValue);
    fAppStarter.CreateTimedApplication(time);
  }
  if (pCmd==fPostRunCmd) {
    fAppStarter.PostRun();
  }
  if (pCmd==fRunCmd) {
    G4int nevents = fRunCmd->GetNewIntValue(szValue);
    fAppStarter.Run(nevents);
  }
}
