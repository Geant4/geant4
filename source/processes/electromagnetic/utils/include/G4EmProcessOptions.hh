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
// $Id: G4EmProcessOptions.hh 96088 2016-03-14 16:03:38Z gcosmo $
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4EmProcessOptions
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 27.02.2004
//
// Modifications:
// 10-01-06 PreciseRange -> CSDARange (V.Ivantchenko)
// 22-05-06 Add SetBremsstrahlungTh (V.Ivanchenko)
// 12-02-07 Add SetSkin, SetLinearLossLimit (V.Ivanchenko)
//
//
// Class Description:
//
// Provide options for EM processes

// -------------------------------------------------------------------
//

#ifndef G4EmProcessOptions_h
#define G4EmProcessOptions_h 1

#include <vector>
#include "globals.hh"
#include "G4MscStepLimitType.hh"

class G4EmParameters;
class G4LossTableManager;

class G4EmProcessOptions
{

public:

  G4EmProcessOptions();

  ~G4EmProcessOptions();

  void SetLossFluctuations(G4bool val);

  void SetBuildCSDARange(G4bool val);

  void SetLPMFlag(G4bool val);

  void SetSplineFlag(G4bool val);

  void SetUseCutAsFinalRange(G4bool val);

  void SetApplyCuts(G4bool val);

  void SetFluo(G4bool val);

  void SetAuger(G4bool val);

  void SetPIXE(G4bool val);

  void SetDeexcitationIgnoreCuts(G4bool val);

  void SetMscLateralDisplacement(G4bool val);

  void SetMscMuHadLateralDisplacement(G4bool val);

  void SetDisplacementBeyondSafety(G4bool val);

  void SetMinSubRange(G4double val);

  void SetMinEnergy(G4double val);

  void SetMaxEnergy(G4double val);

  void SetMaxEnergyForMuons(G4double val);

  void SetMaxEnergyForCSDARange(G4double val);

  void SetLinearLossLimit(G4double val);

  void SetBremsstrahlungTh(G4double val);

  void SetLambdaFactor(G4double val);

  void SetFactorForAngleLimit(G4double val);

  void SetPolarAngleLimit(G4double val);

  void SetMscRangeFactor(G4double val);

  void SetMscGeomFactor(G4double val);

  void SetSkin(G4double val);

  void SetDEDXBinning(G4int val);

  void SetDEDXBinningForCSDARange(G4int val);

  void SetLambdaBinning(G4int val);

  void SetVerbose(G4int val);

  void SetWorkerVerbose(G4int val);

  void SetMscStepLimitation(G4MscStepLimitType val);

  void SetSubCutoff(G4bool val, const G4String& r = "");

  void SetIntegral(G4bool val);

  void SetStepFunction(G4double v1, G4double v2);

  void SetDeexcitationActiveRegion(const G4String& rname = "", 
				   G4bool valDeexcitation = true,
				   G4bool valAuger = true,
				   G4bool valPIXE = true);

  void SetPIXECrossSectionModel(const G4String& val);

  void SetPIXEElectronCrossSectionModel(const G4String& val);

  void SetProcessBiasingFactor(const G4String& name, G4double val, 
			       G4bool flag = true);

  void ActivateForcedInteraction(const G4String& name, G4double length=0.0, 
				 const G4String& region="",
				 G4bool flag = true);

  void ActivateSecondaryBiasing(const G4String& name, const G4String& region, 
				G4double factor, G4double energyLimit);

  void ActivateSecondaryBiasingForGamma(const G4String& name, 
					const G4String& region, 
					G4double factor, G4double energyLimit);

private:

  G4EmProcessOptions & operator=(const  G4EmProcessOptions &right) = delete;
  G4EmProcessOptions(const  G4EmProcessOptions&) = delete;

  G4EmParameters* theParameters;
};

//....oooOO0OOooo.......oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
