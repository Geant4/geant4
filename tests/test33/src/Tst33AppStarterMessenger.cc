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
// $Id: Tst33AppStarterMessenger.cc,v 1.3 2002-10-31 08:32:44 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33AppStarterMessenger.cc
//
// ----------------------------------------------------------------------

#include "Tst33AppStarterMessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "Tst33AppStarter.hh"

Tst33AppStarterMessenger::
Tst33AppStarterMessenger(Tst33AppStarter &appstarter)
  :
  G4UImessenger(),
  fAppStarter(appstarter),
  fMassGeoCmd(new G4UIcommand("/Tst33/MassGeometry", this)),
  fParallelGeoCmd(new G4UIcommand("/Tst33/ParallelGeometry",this)),
  fScoringCmd(new G4UIcommand("/Tst33/Scoring",this)),
  fImpCmd(new G4UIcommand("/Tst33/ImportanceSampling",this)),
  fWWRCmd(new G4UIcommand("/Tst33/WeightRoulette",this)),
  fClearSmaplingCmd(new G4UIcommand("/Tst33/ClearSampling",this)),
  fConfigureSamplingCmd(new G4UIcommand("/Tst33/ConfigureSampling",this)),
  fVisAppComand(new G4UIcommand("/Tst33/Visualization",this)),
  fTimedAppComand(new G4UIcmdWithAnInteger("/Tst33/Timed",this)),
  fPostRunCmd(new G4UIcommand("/Tst33/PostRun",this)),
  fRunCmd(new G4UIcmdWithAnInteger("/Tst33/Run", this))
{
}

Tst33AppStarterMessenger::
~Tst33AppStarterMessenger(){
  delete fMassGeoCmd;
  delete fParallelGeoCmd;
  delete fScoringCmd;
  delete fImpCmd;
  delete fWWRCmd;
  delete fClearSmaplingCmd;
  delete fConfigureSamplingCmd;
  delete fVisAppComand;
  delete fTimedAppComand;
  delete fPostRunCmd;
  delete fRunCmd;
}


void Tst33AppStarterMessenger::SetNewValue(G4UIcommand* pCmd,
					 G4String szValue) {
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
  if (pCmd==fWWRCmd) {
    fAppStarter.CreateWeightRoulette();
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
    G4int time(fTimedAppComand->GetNewIntValue(szValue));
    fAppStarter.CreateTimedApplication(time);
  }
  if (pCmd==fPostRunCmd) {
    fAppStarter.PostRun();
  }
  if (pCmd==fRunCmd) {
    G4int nevents(fRunCmd->GetNewIntValue(szValue));
    fAppStarter.Run(nevents);
  }
  
}



