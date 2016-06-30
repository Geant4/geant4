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
// $Id: G4EmProcessOptions.cc 96088 2016-03-14 16:03:38Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4EmProcessOptions
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 27.02.2004
//
// Modifications:
// 30-06-04 G4EmProcess is pure discrete (V.Ivanchenko)
// 24-03-05 Add ApplyCuts and RandomStep (V.Ivanchenko)
// 10-01-06 PreciseRange -> CSDARange (V.Ivantchenko)
// 10-05-06 Add command MscStepLimit to G4LossTableManager (V.Ivantchenko) 
// 22-05-06 Add SetBremsstrahlungTh (V.Ivanchenko)
// 12-02-07 Add SetSkin, SetLinearLossLimit (V.Ivanchenko)
// 30-05-12 Add biasing for G4VEmProcess (D. Sawkey)
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmProcessOptions.hh"
#include "G4EmParameters.hh"
#include "G4SystemOfUnits.hh"
#include "G4VEmProcess.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4Threading.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmProcessOptions::G4EmProcessOptions()
{
  G4cout << "### WARNING: G4EmProcessOptions class is obsolete and "
	 << "will be removed in the next public release \n"
	 << "    Please, try to use G4EmParameters class and/or UI "
	 << "interface to EM parameters" 
	 << G4endl;
  theParameters = G4EmParameters::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmProcessOptions::~G4EmProcessOptions()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetLossFluctuations(G4bool val)
{
  theParameters->SetLossFluctuations(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetBuildCSDARange(G4bool val)
{
  theParameters->SetBuildCSDARange(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetLPMFlag(G4bool val)
{
  theParameters->SetLPM(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetSplineFlag(G4bool val)
{
  theParameters->SetSpline(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetUseCutAsFinalRange(G4bool val)
{
  theParameters->SetUseCutAsFinalRange(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetApplyCuts(G4bool val)
{
  theParameters->SetApplyCuts(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetFluo(G4bool val)
{
  theParameters->SetFluo(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetAuger(G4bool val)
{
  theParameters->SetAuger(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetPIXE(G4bool val)
{
  theParameters->SetPixe(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetDeexcitationIgnoreCuts(G4bool val)
{
  theParameters->SetDeexcitationIgnoreCut(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMscLateralDisplacement(G4bool val)
{
  theParameters->SetLateralDisplacement(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMscMuHadLateralDisplacement(G4bool val)
{
  theParameters->SetMuHadLateralDisplacement(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetDisplacementBeyondSafety(G4bool val)
{
  theParameters->SetLatDisplacementBeyondSafety(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMinSubRange(G4double val)
{
  theParameters->SetMinSubRange(val); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMinEnergy(G4double val)
{
  theParameters->SetMinEnergy(val); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMaxEnergy(G4double val)
{
  theParameters->SetMaxEnergy(val); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMaxEnergyForMuons(G4double val)
{
  theParameters->SetMaxEnergy(val); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMaxEnergyForCSDARange(G4double val)
{
  theParameters->SetMaxEnergyForCSDARange(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetLinearLossLimit(G4double val)
{
  theParameters->SetLinearLossLimit(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetBremsstrahlungTh(G4double val)
{
  theParameters->SetBremsstrahlungTh(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetLambdaFactor(G4double val)
{
  theParameters->SetLambdaFactor(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetFactorForAngleLimit(G4double val)
{
  theParameters->SetFactorForAngleLimit(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetPolarAngleLimit(G4double val)
{
  theParameters->SetMscThetaLimit(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMscRangeFactor(G4double val)
{
  theParameters->SetMscRangeFactor(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMscGeomFactor(G4double val)
{
  theParameters->SetMscGeomFactor(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetSkin(G4double val)
{
  theParameters->SetMscSkin(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetDEDXBinning(G4int val)
{
  theParameters->SetNumberOfBins(val); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetDEDXBinningForCSDARange(G4int)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetLambdaBinning(G4int val)
{
  theParameters->SetNumberOfBins(val); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetVerbose(G4int val)
{
  theParameters->SetVerbose(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetWorkerVerbose(G4int val)
{
  theParameters->SetWorkerVerbose(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetMscStepLimitation(G4MscStepLimitType val)
{
  theParameters->SetMscStepLimitType(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetSubCutoff(G4bool val, const G4String& r)
{
  theParameters->SetSubCutoff(val, r);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetIntegral(G4bool val)
{
  theParameters->SetIntegral(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetStepFunction(G4double v1, G4double v2)
{
  theParameters->SetStepFunction(v1, v2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4EmProcessOptions::SetDeexcitationActiveRegion(const G4String& rname, 
						G4bool valDeex,
						G4bool valAuger,
						G4bool valPIXE)
{
  theParameters->SetDeexActiveRegion(rname, valDeex, valAuger, valPIXE);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmProcessOptions::SetPIXECrossSectionModel(const G4String& mname)
{
  theParameters->SetPIXECrossSectionModel(mname); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4EmProcessOptions::SetPIXEElectronCrossSectionModel(const G4String& mname)
{
  theParameters->SetPIXEElectronCrossSectionModel(mname); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4EmProcessOptions::SetProcessBiasingFactor(const G4String& name, G4double val,
					    G4bool flag)
{
  theParameters->SetProcessBiasingFactor(name, val, flag);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4EmProcessOptions::ActivateForcedInteraction(const G4String& name, 
					      G4double length, 
					      const G4String& region,
					      G4bool flag)
{
  theParameters->ActivateForcedInteraction(name, region, length, flag);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4EmProcessOptions::ActivateSecondaryBiasing(const G4String& name,
					     const G4String& region, 
					     G4double factor,
					     G4double energyLimit)
{
  theParameters->ActivateSecondaryBiasing(name, region, factor, energyLimit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4EmProcessOptions::ActivateSecondaryBiasingForGamma(const G4String& name,
					     const G4String& region, 
					     G4double factor,
					     G4double energyLimit)
{
  theParameters->ActivateSecondaryBiasing(name, region, factor, energyLimit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

