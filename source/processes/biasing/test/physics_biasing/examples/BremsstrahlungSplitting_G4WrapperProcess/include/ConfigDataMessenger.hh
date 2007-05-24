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
// $Id: ConfigDataMessenger.hh,v 1.1 2007-05-24 21:57:03 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation - messenger for ConfigData
//
#ifndef CONFIGDATAMESSENGER_HH
#define CONFIGDATAMESSENGER_HH 

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIdirectory;

class ConfigDataMessenger : public G4UImessenger {

public:
  
  ConfigDataMessenger();
  
  virtual ~ConfigDataMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  
  // Built in configurations
  G4UIcmdWithoutParameter* fBe_15pt8MeV_1_10_degrees;
  G4UIcmdWithoutParameter* fBe_15pt18MeV_30_90_degrees;

  G4UIcmdWithoutParameter* fAl_15pt8MeV_1_10_degrees;
  G4UIcmdWithoutParameter* fAl_15pt18MeV_30_90_degrees;

  G4UIcmdWithoutParameter* fPb_15pt8MeV_1_10_degrees;
  G4UIcmdWithoutParameter* fPb_15pt18MeV_30_90_degrees;

  ////////////////////
  G4UIcmdWithoutParameter* fPb_10pt09MeV_0_degrees;
  G4UIcmdWithoutParameter* fPb_15pt18MeV_0_degrees;
  G4UIcmdWithoutParameter* fPb_20pt28MeV_0_degrees;
  G4UIcmdWithoutParameter* fPb_25pt38MeV_0_degrees;
  G4UIcmdWithoutParameter* fPb_30pt45MeV_0_degrees;

  G4UIcmdWithoutParameter* fAl_10pt09MeV_0_degrees;
  G4UIcmdWithoutParameter* fAl_15pt18MeV_0_degrees;
  G4UIcmdWithoutParameter* fAl_20pt28MeV_0_degrees;
  G4UIcmdWithoutParameter* fAl_25pt38MeV_0_degrees;
  G4UIcmdWithoutParameter* fAl_30pt45MeV_0_degrees;

  G4UIcmdWithoutParameter* fVerbose;

  G4UIdirectory* listDir;

  G4UIcmdWithADoubleAndUnit* scorerMinRadiusCmd;
  G4UIcmdWithADoubleAndUnit* scorerMaxRadiusCmd;
  G4UIcmdWithADoubleAndUnit* scorerMinThetaCmd;
  G4UIcmdWithADoubleAndUnit* scorerMaxThetaCmd;
  G4UIcmdWithADoubleAndUnit* scorerDeltaThetaCmd;
  G4UIcmdWithADoubleAndUnit* scorerMinEnergyCmd;
  G4UIcmdWithADoubleAndUnit* scorerMaxEnergyCmd;
  G4UIcmdWithADoubleAndUnit* scorerDeltaEnergyCmd;

  G4UIcmdWithAnInteger* bremSplittingFactor;
  G4UIcmdWithABool* bremSplittingActivation;

};

#endif








