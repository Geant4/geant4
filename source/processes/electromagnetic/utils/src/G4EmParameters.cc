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
// $Id: G4EmParameters.cc 69320 2013-04-30 15:59:36Z vnivanch $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4EmParameters
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 18.05.2013
//
// Modifications:
//
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmParameters.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmParametersMessenger.hh"
#include "G4NistManager.hh"

G4EmParameters* G4EmParameters::theInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmParameters* G4EmParameters::Instance()
{
  if(0 == theInstance) {
    static G4EmParameters manager;
    theInstance = &manager;
  }
  return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmParameters::~G4EmParameters()
{
  delete theMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmParameters::G4EmParameters()
{
  G4NistManager::Instance();
  theMessenger = new G4EmParametersMessenger(this);

  lossFluctuation = true;
  buildCSDARange = false;
  flagLPM = true;
  spline = true;
  finalRange = false;
  applyCuts = false;
  fluo = false;
  auger = false;
  pixe = false;
  deexIgnoreCut = false;
  lateralDisplacement = true;
  muhadLateralDisplacement = true;
  latDisplacementBeyondSafety = false;
  useAngGeneratorForIonisation = false;

  minSubRange = DBL_MAX;
  minKinEnergy = 0.1*keV;
  maxKinEnergy = 10.0*TeV;
  maxKinEnergyCSDA = 1.0*GeV;
  linLossLimit = 0.01;
  bremsTh = DBL_MAX;
  lambdaFactor = 0.8;
  factorForAngleLimit = 1.0;
  thetaLimit = CLHEP::pi;
  rangeFactor = 0.04;
  geomFactor = 2.5;
  skin = 1.0;

  nbins  = 77;
  nbinsPerDecade = 7;
  verbose = 1;
  workerVerbose = 0;

  mscStepLimit = fUseSafety;
}

#include "G4AutoLock.hh"
namespace { G4Mutex EmParametersMutex = G4MUTEX_INITIALIZER; }

void G4EmParameters::SetLossFluctuations(G4bool val)
{
  G4AutoLock l(&EmParametersMutex);
  lossFluctuation = val;
}

G4bool G4EmParameters::LossFluctuation() const
{
  return lossFluctuation;
}

void G4EmParameters::SetBuildCSDARange(G4bool val)
{
  G4AutoLock l(&EmParametersMutex);
  buildCSDARange = val;
}

G4bool G4EmParameters::BuildCSDARange() const 
{
  return buildCSDARange;
}

void G4EmParameters::SetLPM(G4bool val)
{
  G4AutoLock l(&EmParametersMutex);
  flagLPM = val;
}

G4bool G4EmParameters::LPM() const 
{
  return flagLPM;
}

void G4EmParameters::SetSpline(G4bool val)
{
  G4AutoLock l(&EmParametersMutex);
  spline = val;
}

G4bool G4EmParameters::Spline() const
{
  return spline;
}

void G4EmParameters::SetUseCutAsFinalRange(G4bool val)
{
  G4AutoLock l(&EmParametersMutex);
  finalRange = val;
}

G4bool G4EmParameters::UseCutAsFinalRange() const
{
  return finalRange;
}

void G4EmParameters::SetApplyCuts(G4bool val)
{
  G4AutoLock l(&EmParametersMutex);
  applyCuts = val;
}

G4bool G4EmParameters::ApplyCuts() const
{
  return applyCuts;
}

void G4EmParameters::SetFluo(G4bool val)
{
  G4AutoLock l(&EmParametersMutex);
  fluo = val;
}

G4bool G4EmParameters::Fluo() const
{
  return fluo;
}

void G4EmParameters::SetAuger(G4bool val)
{
  G4AutoLock l(&EmParametersMutex);
  auger = val;
  if(val) { fluo = true; }
}

G4bool G4EmParameters::Auger() const
{
  return auger;
}

void G4EmParameters::SetPixe(G4bool val)
{
  G4AutoLock l(&EmParametersMutex);
  pixe = val;
  if(val) { fluo = true; }
}

G4bool G4EmParameters::Pixe() const
{
  return pixe;
}

void G4EmParameters::SetDeexcitationIgnoreCut(G4bool val)
{
  G4AutoLock l(&EmParametersMutex);
  deexIgnoreCut = val;
}

G4bool G4EmParameters::DeexcitationIgnoreCut() const
{
  return deexIgnoreCut;
}

void G4EmParameters::SetLateralDisplacement(G4bool val)
{
  G4AutoLock l(&EmParametersMutex);
  lateralDisplacement = val;
}

G4bool G4EmParameters::LateralDisplacement() const
{
  return lateralDisplacement;
}

void G4EmParameters::SetMuHadLateralDisplacement(G4bool val)
{
  G4AutoLock l(&EmParametersMutex);
  muhadLateralDisplacement = val;
}

G4bool G4EmParameters::MuHadLateralDisplacement() const
{
  return muhadLateralDisplacement;
}

void G4EmParameters::SetLatDisplacementBeyondSafety(G4bool val)
{
  G4AutoLock l(&EmParametersMutex);
  latDisplacementBeyondSafety = val;
}

G4bool G4EmParameters::LatDisplacementBeyondSafety() const
{
  return latDisplacementBeyondSafety;
}

void G4EmParameters::ActivateAngularGeneratorForIonisation(G4bool val)
{
  useAngGeneratorForIonisation = val;
}

G4bool G4EmParameters::UseAngularGeneratorForIonisation() const
{
  return useAngGeneratorForIonisation;
}

void G4EmParameters::SetMinSubRange(G4double val)
{
  G4AutoLock l(&EmParametersMutex);
  if(val > 0.0 && val < 1.0) {
    minSubRange = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of MinSubRange is out of range (0 - 1): " << val
       << " is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::MinSubRange() const
{
  return minSubRange;
}

void G4EmParameters::SetMinEnergy(G4double val)
{
  G4AutoLock l(&EmParametersMutex);
  if(val > 1.e-3*eV && val < maxKinEnergy) {
    minKinEnergy = val;
    nbins = nbinsPerDecade*G4lrint(std::log10(maxKinEnergy/minKinEnergy));
  } else {
    G4ExceptionDescription ed;
    ed << "Value of MinKinEnergy is out of range: " << val/MeV 
       << " MeV is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::MinKinEnergy() const
{
  return minKinEnergy;
}

void G4EmParameters::SetMaxEnergy(G4double val)
{
  G4AutoLock l(&EmParametersMutex);
  if(val > minKinEnergy && val < 1.e+7*TeV) {
    maxKinEnergy = val;
    nbins = nbinsPerDecade*G4lrint(std::log10(maxKinEnergy/minKinEnergy));
  } else {
    G4ExceptionDescription ed;
    ed << "Value of MaxKinEnergy is out of range: " 
       << val/GeV << " GeV is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::MaxKinEnergy() const
{
  return maxKinEnergy;
}

void G4EmParameters::SetMaxEnergyForCSDARange(G4double val)
{
  G4AutoLock l(&EmParametersMutex);
  if(val > minKinEnergy && val <= 100*TeV) {
    maxKinEnergyCSDA = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of MaxKinEnergyCSDA is out of range: " 
       << val/GeV << " GeV is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::MaxEnergyForCSDARange() const
{
  return maxKinEnergyCSDA; 
}

void G4EmParameters::SetLinearLossLimit(G4double val)
{
  G4AutoLock l(&EmParametersMutex);
  if(val > 0.0 && val < 0.5) {
    linLossLimit = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of linLossLimit is out of range: " << val 
       << " is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::LinearLossLimit() const
{
  return linLossLimit;
}

void G4EmParameters::SetBremsstrahlungTh(G4double val)
{
  G4AutoLock l(&EmParametersMutex);
  if(val > 0.0) {
    bremsTh = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of bremsstrahlung threshold is out of range: " 
       << val/GeV << " GeV is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::BremsstrahlungTh() const 
{
  return bremsTh;
}

void G4EmParameters::SetLambdaFactor(G4double val)
{
  G4AutoLock l(&EmParametersMutex);
  if(val > 0.0 && val < 1.0) {
    lambdaFactor = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of lambda factor is out of range: " << val 
       << " is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::LambdaFactor() const 
{
  return lambdaFactor;
}

void G4EmParameters::SetFactorForAngleLimit(G4double val)
{
  G4AutoLock l(&EmParametersMutex);
  if(val > 0.0) {
    factorForAngleLimit = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of factor for enegry limit is out of range: " 
       << val << " is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::FactorForAngleLimit() const 
{
  return factorForAngleLimit;
}

void G4EmParameters::SetMscThetaLimit(G4double val)
{
  G4AutoLock l(&EmParametersMutex);
  if(val >= 0.0 && val <= pi) {
    thetaLimit = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of polar angle limit is out of range: " 
       << val << " is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::MscThetaLimit() const 
{
  return thetaLimit;
}

void G4EmParameters::SetMscRangeFactor(G4double val)
{
  G4AutoLock l(&EmParametersMutex);
  if(val > 0.0 && val < 1.0) {
    rangeFactor = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of rangeFactor is out of range: " 
       << val << " is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::MscRangeFactor() const 
{
  return rangeFactor;
}

void G4EmParameters::SetMscGeomFactor(G4double val)
{
  G4AutoLock l(&EmParametersMutex);
  if(val >= 1.0) {
    geomFactor = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of geomFactor is out of range: " 
       << val << " is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::MscGeomFactor() const 
{
  return geomFactor;
}

void G4EmParameters::SetMscSkin(G4double val)
{
  G4AutoLock l(&EmParametersMutex);
  if(val >= 1.0) {
    skin = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of skin is out of range: " 
       << val << " is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::MscSkin() const 
{
  return skin;
}

void G4EmParameters::SetNumberOfBins(G4int val)
{
  G4AutoLock l(&EmParametersMutex);
  if(val >= 5 && val < 10000000) {
    nbins = val;
    nbinsPerDecade = G4lrint(nbins/std::log10(maxKinEnergy/minKinEnergy));
  } else {
    G4ExceptionDescription ed;
    ed << "Value of number of bins is out of range: " 
       << val << " is ignored"; 
    PrintWarning(ed);
  }
}

G4int G4EmParameters::NumberOfBins() const 
{
  return nbins;
}

void G4EmParameters::SetNumberOfBinsPerDecade(G4int val)
{
  G4AutoLock l(&EmParametersMutex);
  if(val >= 5 && val < 1000000) {
    nbinsPerDecade = val;
    nbins = nbinsPerDecade*G4lrint(std::log10(maxKinEnergy/minKinEnergy));
  } else {
    G4ExceptionDescription ed;
    ed << "Value of number of bins per decade is out of range: " 
       << val << " is ignored"; 
    PrintWarning(ed);
  }
}

G4int G4EmParameters::NumberOfBinsPerDecade() const 
{
  return nbinsPerDecade; 
}

void G4EmParameters::SetVerbose(G4int val)
{
  G4AutoLock l(&EmParametersMutex);
  verbose = val;
  workerVerbose = std::min(workerVerbose, verbose);
}

G4int G4EmParameters::Verbose() const 
{
  return verbose;
}

void G4EmParameters::SetWorkerVerbose(G4int val)
{
  G4AutoLock l(&EmParametersMutex);
  workerVerbose = val;
}

G4int G4EmParameters::WorkerVerbose() const 
{
  return workerVerbose;
}

void G4EmParameters::SetMscStepLimitType(G4MscStepLimitType val)
{
  G4AutoLock l(&EmParametersMutex);
  mscStepLimit = val;
}

G4MscStepLimitType G4EmParameters::MscStepLimitType() const 
{
  return mscStepLimit;
}

void G4EmParameters::PrintWarning(G4ExceptionDescription& ed)
{
  G4Exception("G4EmParameters", "em0044", JustWarning, ed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

