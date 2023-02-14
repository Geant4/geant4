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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4EmParameters
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 18.05.2013
//
// Modifications:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmParameters.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VEmProcess.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4EmExtraParameters.hh"
#include "G4EmLowEParameters.hh"
#include "G4EmParametersMessenger.hh"
#include "G4NistManager.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ApplicationState.hh"
#include "G4StateManager.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

G4EmParameters* G4EmParameters::theInstance = nullptr;

namespace
{
  G4Mutex emParametersMutex = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmParameters* G4EmParameters::Instance()
{
  if(nullptr == theInstance) { 
    G4AutoLock l(&emParametersMutex);
    if(nullptr == theInstance) {
      static G4EmParameters manager;
      theInstance = &manager;
    }
    l.unlock();
  }
  return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmParameters::~G4EmParameters()
{
  delete theMessenger;
  delete fBParameters;
  delete fCParameters;
  delete emSaturation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmParameters::G4EmParameters()
{
  G4NistManager::Instance();
  theMessenger = new G4EmParametersMessenger(this);
  Initialise();

  fBParameters = new G4EmExtraParameters();
  fCParameters = new G4EmLowEParameters();

  fStateManager = G4StateManager::GetStateManager();
  emSaturation = nullptr;
}

void G4EmParameters::SetDefaults()
{
  if(!IsLocked()) { 
    Initialise();
    fBParameters->Initialise();
    fCParameters->Initialise();
  }
}

void G4EmParameters::Initialise()
{
  lossFluctuation = true;
  buildCSDARange = false;
  flagLPM = true;
  cutAsFinalRange = false;
  applyCuts = false;
  lateralDisplacement = true;
  lateralDisplacementAlg96 = true;
  muhadLateralDisplacement = false;
  useAngGeneratorForIonisation = false;
  useMottCorrection = false;
  integral = true;
  birks = false;
  fICRU90 = false;
  gener = false;
  onIsolated = false;
  fSamplingTable = false;
  fPolarisation = false;
  fMuDataFromFile = false;
  fPEKShell = true;
  fMscPosiCorr = true;
  fDNA = false;
  fIsPrinted = false;

  minKinEnergy = 0.1*CLHEP::keV;
  maxKinEnergy = 100.0*CLHEP::TeV;
  maxKinEnergyCSDA = 1.0*CLHEP::GeV;
  max5DEnergyForMuPair = 0.0;
  lowestElectronEnergy = 1.0*CLHEP::keV;
  lowestMuHadEnergy = 1.0*CLHEP::keV;
  lowestTripletEnergy = 1.0*CLHEP::MeV;
  maxNIELEnergy = 0.0;
  linLossLimit = 0.01;
  bremsTh = bremsMuHadTh = maxKinEnergy;
  lambdaFactor = 0.8;
  factorForAngleLimit = 1.0;
  thetaLimit = CLHEP::pi;
  energyLimit = 100.0*CLHEP::MeV;
  rangeFactor = 0.04;
  rangeFactorMuHad = 0.2;
  geomFactor = 2.5;
  skin = 1.0;
  safetyFactor = 0.6;
  lambdaLimit  = 1.0*CLHEP::mm;
  factorScreen = 1.0;

  nbinsPerDecade = 7;
  verbose = 1;
  workerVerbose = 0;
  tripletConv = 0;

  fTransportationWithMsc = G4TransportationWithMscType::fDisabled;
  mscStepLimit = fUseSafety;
  mscStepLimitMuHad = fMinimal;
  nucFormfactor = fExponentialNF;
  fSStype = fWVI;
  fFluct = fUniversalFluctuation;
}

void G4EmParameters::SetLossFluctuations(G4bool val)
{
  if(IsLocked()) { return; }
  lossFluctuation = val;
}

G4bool G4EmParameters::LossFluctuation() const
{
  return lossFluctuation;
}

void G4EmParameters::SetBuildCSDARange(G4bool val)
{
  if(IsLocked()) { return; }
  buildCSDARange = val;
}

G4bool G4EmParameters::BuildCSDARange() const 
{
  return buildCSDARange;
}

void G4EmParameters::SetLPM(G4bool val)
{
  if(IsLocked()) { return; }
  flagLPM = val;
}

G4bool G4EmParameters::LPM() const 
{
  return flagLPM;
}

void G4EmParameters::SetUseCutAsFinalRange(G4bool val)
{
  if(IsLocked()) { return; }
  cutAsFinalRange = val;
}

G4bool G4EmParameters::UseCutAsFinalRange() const
{
  return cutAsFinalRange;
}

void G4EmParameters::SetApplyCuts(G4bool val)
{
  if(IsLocked()) { return; }
  applyCuts = val;
}

G4bool G4EmParameters::ApplyCuts() const
{
  return applyCuts;
}

void G4EmParameters::SetFluo(G4bool val)
{
  if(IsLocked()) { return; }
  fCParameters->SetFluo(val);
}

G4bool G4EmParameters::Fluo() const
{
  return fCParameters->Fluo();
}

G4EmFluoDirectory G4EmParameters::FluoDirectory() const
{
  return fCParameters->FluoDirectory();
}

void G4EmParameters::SetFluoDirectory(G4EmFluoDirectory val)
{
  if(IsLocked()) { return; }
  fCParameters->SetFluoDirectory(val);
}

void G4EmParameters::SetBeardenFluoDir(G4bool val)
{
  if(IsLocked()) { return; }
  fCParameters->SetBeardenFluoDir(val);
}

void G4EmParameters::SetANSTOFluoDir(G4bool val)
{
  if(IsLocked()) { return; }
  fCParameters->SetANSTOFluoDir(val);
}

void G4EmParameters::SetXDB_EADLFluoDir(G4bool val)
{
  if(IsLocked()) { return; }
  fCParameters->SetXDB_EADLFluoDir(val);
}

void G4EmParameters::SetAuger(G4bool val)
{
  if(IsLocked()) { return; }
  fCParameters->SetAuger(val);
}

G4bool G4EmParameters::BeardenFluoDir()
{
  auto dir = fCParameters->FluoDirectory();
  return (dir == fluoBearden);
}

G4bool G4EmParameters::ANSTOFluoDir()
{
  auto dir = fCParameters->FluoDirectory();
  return (dir == fluoANSTO);
}

G4bool G4EmParameters::Auger() const
{
  return fCParameters->Auger();
}

void G4EmParameters::SetPixe(G4bool val)
{
  if(IsLocked()) { return; }
  fCParameters->SetPixe(val);
}

G4bool G4EmParameters::Pixe() const
{
  return fCParameters->Pixe();
}

void G4EmParameters::SetDeexcitationIgnoreCut(G4bool val)
{
  if(IsLocked()) { return; }
  fCParameters->SetDeexcitationIgnoreCut(val);
}

G4bool G4EmParameters::DeexcitationIgnoreCut() const
{
  return fCParameters->DeexcitationIgnoreCut();
}

void G4EmParameters::SetLateralDisplacement(G4bool val)
{
  if(IsLocked()) { return; }
  lateralDisplacement = val;
}

G4bool G4EmParameters::LateralDisplacement() const
{
  return lateralDisplacement;
}

void G4EmParameters::SetLateralDisplacementAlg96(G4bool val)
{
  if(IsLocked()) { return; }
  lateralDisplacementAlg96 = val;
}

G4bool G4EmParameters::LateralDisplacementAlg96() const
{
  return lateralDisplacementAlg96;
}

void G4EmParameters::SetMuHadLateralDisplacement(G4bool val)
{
  if(IsLocked()) { return; }
  muhadLateralDisplacement = val;
}

G4bool G4EmParameters::MuHadLateralDisplacement() const
{
  return muhadLateralDisplacement;
}

void G4EmParameters::ActivateAngularGeneratorForIonisation(G4bool val)
{
  if(IsLocked()) { return; }
  useAngGeneratorForIonisation = val;
}

G4bool G4EmParameters::UseAngularGeneratorForIonisation() const
{
  return useAngGeneratorForIonisation;
}

void G4EmParameters::SetUseMottCorrection(G4bool val)
{
  if(IsLocked()) { return; }
  useMottCorrection = val;
}

G4bool G4EmParameters::UseMottCorrection() const
{
  return useMottCorrection;
}

void G4EmParameters::SetIntegral(G4bool val)
{
  if(IsLocked()) { return; }
  integral = val;
}

G4bool G4EmParameters::Integral() const
{
  return integral;
}

void G4EmParameters::SetEnablePolarisation(G4bool val)
{
  if(IsLocked()) { return; }
  fPolarisation = val;
}

G4bool G4EmParameters::EnablePolarisation() const
{
  return fPolarisation;
}

void G4EmParameters::SetBirksActive(G4bool val)
{
  if(IsLocked()) { return; }
  birks = val;
  if(birks && nullptr == emSaturation) { emSaturation = new G4EmSaturation(1); }
}

G4bool G4EmParameters::BirksActive() const
{
  return birks;
}

void G4EmParameters::SetUseICRU90Data(G4bool val)
{
  if(IsLocked()) { return; }
  fICRU90 = val;
}

G4bool G4EmParameters::UseICRU90Data() const
{
  return fICRU90;
}

void G4EmParameters::SetDNAFast(G4bool val)
{
  if(IsLocked()) { return; }
  fCParameters->SetDNAFast(val);
  if(val) { ActivateDNA(); }
}

G4bool G4EmParameters::DNAFast() const
{
  return fCParameters->DNAFast();
}

void G4EmParameters::SetDNAStationary(G4bool val)
{
  if(IsLocked()) { return; }
  fCParameters->SetDNAStationary(val);
  if(val) { ActivateDNA(); }
}

G4bool G4EmParameters::DNAStationary() const
{
  return fCParameters->DNAStationary();
}

void G4EmParameters::SetDNAElectronMsc(G4bool val)
{
  if(IsLocked()) { return; }
  fCParameters->SetDNAElectronMsc(val);
  if(val) { ActivateDNA(); }
}

G4bool G4EmParameters::DNAElectronMsc() const
{
  return fCParameters->DNAElectronMsc();
}

void G4EmParameters::SetGeneralProcessActive(G4bool val)
{
  if(IsLocked()) { return; }
  gener = val;
}

G4bool G4EmParameters::GeneralProcessActive() const
{
  return gener;
}

void G4EmParameters::SetEmSaturation(G4EmSaturation* ptr)
{
  if(IsLocked()) { return; }
  birks = (nullptr != ptr);
  if(emSaturation != ptr) {
    delete emSaturation;
    emSaturation = ptr;
  }
}

G4bool G4EmParameters::RetrieveMuDataFromFile() const
{
  return fMuDataFromFile;
}

void G4EmParameters::SetRetrieveMuDataFromFile(G4bool v)
{
  fMuDataFromFile = v;
}

void G4EmParameters::SetOnIsolated(G4bool val)
{
  if(IsLocked()) { return; }
  onIsolated = val;
}

G4bool G4EmParameters::OnIsolated() const
{
  return onIsolated;
}

void G4EmParameters::SetEnableSamplingTable(G4bool val)
{
  if(IsLocked()) { return; }
  fSamplingTable = val;
}

G4bool G4EmParameters::EnableSamplingTable() const
{
  return fSamplingTable;
}

G4bool G4EmParameters::PhotoeffectBelowKShell() const
{
  return fPEKShell;
}

void G4EmParameters::SetPhotoeffectBelowKShell(G4bool v)
{
  if(IsLocked()) { return; }
  fPEKShell = v;
}

G4bool G4EmParameters::MscPositronCorrection() const
{
  return fMscPosiCorr;
}

void G4EmParameters::SetMscPositronCorrection(G4bool v)
{
  if(IsLocked()) { return; }
  fMscPosiCorr = v;
}

void G4EmParameters::ActivateDNA()
{
  if(IsLocked()) { return; }
  fDNA = true;
}

void G4EmParameters::SetIsPrintedFlag(G4bool val)
{
  fIsPrinted = val;
}

G4bool G4EmParameters::IsPrintLocked() const
{
  return fIsPrinted;
}

G4EmSaturation* G4EmParameters::GetEmSaturation()
{
  if(nullptr == emSaturation) { 
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&emParametersMutex);
    if(nullptr == emSaturation) { 
#endif
      emSaturation = new G4EmSaturation(1);
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&emParametersMutex);
#endif
  }
  birks = true;
  return emSaturation;
}

void G4EmParameters::SetMinEnergy(G4double val)
{
  if(IsLocked()) { return; }
  if(val > 1.e-3*CLHEP::eV && val < maxKinEnergy) {
    minKinEnergy = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of MinKinEnergy - is out of range: " << val/CLHEP::MeV 
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
  if(IsLocked()) { return; }
  if(val > std::max(minKinEnergy,9.99*CLHEP::MeV) && val < 1.e+7*CLHEP::TeV) {
    maxKinEnergy = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of MaxKinEnergy is out of range: " 
       << val/CLHEP::GeV 
       << " GeV is ignored; allowed range 10 MeV - 1.e+7 TeV"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::MaxKinEnergy() const
{
  return maxKinEnergy;
}

void G4EmParameters::SetMaxEnergyForCSDARange(G4double val)
{
  if(IsLocked()) { return; }
  if(val > minKinEnergy && val <= 100*CLHEP::TeV) {
    maxKinEnergyCSDA = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of MaxKinEnergyCSDA is out of range: " 
       << val/CLHEP::GeV << " GeV is ignored; allowed range "
       << minKinEnergy << " MeV - 100 TeV"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::MaxEnergyForCSDARange() const
{
  return maxKinEnergyCSDA; 
}

void G4EmParameters::SetLowestElectronEnergy(G4double val)
{
  if(IsLocked()) { return; }
  if(val >= 0.0) { lowestElectronEnergy = val; }
}

G4double G4EmParameters::LowestElectronEnergy() const
{
  return lowestElectronEnergy; 
}

void G4EmParameters::SetLowestMuHadEnergy(G4double val)
{
  if(IsLocked()) { return; }
  if(val >= 0.0) { lowestMuHadEnergy = val; }
}

G4double G4EmParameters::LowestMuHadEnergy() const
{
  return lowestMuHadEnergy; 
}

void G4EmParameters::SetLowestTripletEnergy(G4double val)
{
  if(IsLocked()) { return; }
  if(val > 0.0) { lowestTripletEnergy = val; }
}

G4double G4EmParameters::LowestTripletEnergy() const
{
  return lowestTripletEnergy;
}

void G4EmParameters::SetMaxNIELEnergy(G4double val)
{
  if(IsLocked()) { return; }
  if(val >= 0.0) { maxNIELEnergy = val; }
}

G4double G4EmParameters::MaxNIELEnergy() const
{
  return maxNIELEnergy;
}

void G4EmParameters::SetMaxEnergyFor5DMuPair(G4double val)
{
  if(IsLocked()) { return; }
  if(val > 0.0) { max5DEnergyForMuPair = val; }
}

G4double G4EmParameters::MaxEnergyFor5DMuPair() const
{
  return max5DEnergyForMuPair;
}

void G4EmParameters::SetLinearLossLimit(G4double val)
{
  if(IsLocked()) { return; }
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
  if(IsLocked()) { return; }
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

void G4EmParameters::SetMuHadBremsstrahlungTh(G4double val)
{
  if(IsLocked()) { return; }
  if(val > 0.0) {
    bremsMuHadTh = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of bremsstrahlung threshold is out of range: " 
       << val/GeV << " GeV is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::MuHadBremsstrahlungTh() const
{
  return bremsMuHadTh;
}

void G4EmParameters::SetLambdaFactor(G4double val)
{
  if(IsLocked()) { return; }
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
  if(IsLocked()) { return; }
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
  if(IsLocked()) { return; }
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

void G4EmParameters::SetMscEnergyLimit(G4double val)
{
  if(IsLocked()) { return; }
  if(val >= 0.0) {
    energyLimit = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of msc energy limit is out of range: " 
       << val << " is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::MscEnergyLimit() const
{
  return energyLimit;
}

void G4EmParameters::SetMscRangeFactor(G4double val)
{
  if(IsLocked()) { return; }
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

void G4EmParameters::SetMscMuHadRangeFactor(G4double val)
{
  if(IsLocked()) { return; }
  if(val > 0.0 && val < 1.0) {
    rangeFactorMuHad = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of rangeFactorMuHad is out of range: " 
       << val << " is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::MscMuHadRangeFactor() const 
{
  return rangeFactorMuHad;
}

void G4EmParameters::SetMscGeomFactor(G4double val)
{
  if(IsLocked()) { return; }
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

void G4EmParameters::SetMscSafetyFactor(G4double val)
{
  if(IsLocked()) { return; }
  if(val >= 0.1) {
    safetyFactor = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of safetyFactor is out of range: " 
       << val << " is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::MscSafetyFactor() const 
{
  return safetyFactor;
}

void G4EmParameters::SetMscLambdaLimit(G4double val)
{
  if(IsLocked()) { return; }
  if(val >= 0.0) {
    lambdaLimit = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of lambdaLimit is out of range: " 
       << val << " is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::MscLambdaLimit() const 
{
  return lambdaLimit;
}

void G4EmParameters::SetMscSkin(G4double val)
{
  if(IsLocked()) { return; }
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

void G4EmParameters::SetScreeningFactor(G4double val)
{
  if(IsLocked()) { return; }
  if(val > 0.0) {
    factorScreen = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of factorScreen is out of range: " 
       << val << " is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::ScreeningFactor() const
{
  return factorScreen;
}

void G4EmParameters::SetStepFunction(G4double v1, G4double v2)
{
  if(IsLocked()) { return; }
  fBParameters->SetStepFunction(v1, v2);
}

void G4EmParameters::SetStepFunctionMuHad(G4double v1, G4double v2)
{
  if(IsLocked()) { return; }
  fBParameters->SetStepFunctionMuHad(v1, v2);
}

void G4EmParameters::SetStepFunctionLightIons(G4double v1, G4double v2)
{
  if(IsLocked()) { return; }
  fBParameters->SetStepFunctionLightIons(v1, v2);
}

void G4EmParameters::SetStepFunctionIons(G4double v1, G4double v2)
{
  if(IsLocked()) { return; }
  fBParameters->SetStepFunctionIons(v1, v2);
}

void G4EmParameters::FillStepFunction(const G4ParticleDefinition* part, G4VEnergyLossProcess* proc) const
{
  fBParameters->FillStepFunction(part, proc);
}

G4int G4EmParameters::NumberOfBins() const 
{
  return nbinsPerDecade*G4lrint(std::log10(maxKinEnergy/minKinEnergy));
}

void G4EmParameters::SetNumberOfBinsPerDecade(G4int val)
{
  if(IsLocked()) { return; }
  if(val >= 5 && val < 1000000) {
    nbinsPerDecade = val;
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
  if(IsLocked()) { return; }
  verbose = val;
  workerVerbose = std::min(workerVerbose, verbose);
}

G4int G4EmParameters::Verbose() const 
{
  return verbose;
}

void G4EmParameters::SetWorkerVerbose(G4int val)
{
  if(IsLocked()) { return; }
  workerVerbose = val;
}

G4int G4EmParameters::WorkerVerbose() const 
{
  return workerVerbose;
}

void G4EmParameters::SetTransportationWithMsc(G4TransportationWithMscType val)
{
  if(IsLocked()) { return; }
  fTransportationWithMsc = val;
}

G4TransportationWithMscType G4EmParameters::TransportationWithMsc() const
{
  return fTransportationWithMsc;
}

void G4EmParameters::SetFluctuationType(G4EmFluctuationType val)
{
  if(IsLocked()) { return; }
  fFluct = val;
}

G4EmFluctuationType G4EmParameters::FluctuationType() const
{
  return fFluct;
}

void G4EmParameters::SetMscStepLimitType(G4MscStepLimitType val)
{
  if(IsLocked()) { return; }
  mscStepLimit = val;
}

G4MscStepLimitType G4EmParameters::MscStepLimitType() const 
{
  return mscStepLimit;
}

void G4EmParameters::SetMscMuHadStepLimitType(G4MscStepLimitType val)
{
  if(IsLocked()) { return; }
  mscStepLimitMuHad = val;
}

G4MscStepLimitType G4EmParameters::MscMuHadStepLimitType() const 
{
  return mscStepLimitMuHad;
}

void G4EmParameters::SetSingleScatteringType(G4eSingleScatteringType val)
{
  if(IsLocked()) { return; }
  fSStype = val;
}

G4eSingleScatteringType G4EmParameters::SingleScatteringType() const
{
  return fSStype;
}

void 
G4EmParameters::SetNuclearFormfactorType(G4NuclearFormfactorType val)
{
  if(IsLocked()) { return; }
  nucFormfactor = val;
}

G4NuclearFormfactorType G4EmParameters::NuclearFormfactorType() const
{
  return nucFormfactor;
}

void G4EmParameters::SetDNAeSolvationSubType(G4DNAModelSubType val)
{
  if(IsLocked()) { return; }
  fCParameters->SetDNAeSolvationSubType(val);
  ActivateDNA();
}

G4DNAModelSubType G4EmParameters::DNAeSolvationSubType() const
{
  return fCParameters->DNAeSolvationSubType();
}

void G4EmParameters::SetConversionType(G4int val)
{
  if(IsLocked()) { return; }
  tripletConv = val;
}

G4int G4EmParameters::GetConversionType() const
{
  return tripletConv;
}

void G4EmParameters::SetPIXECrossSectionModel(const G4String& sss)
{
  if(IsLocked()) { return; }
  fCParameters->SetPIXECrossSectionModel(sss);
}

const G4String& G4EmParameters::PIXECrossSectionModel()
{
  return fCParameters->PIXECrossSectionModel();
}

void G4EmParameters::SetPIXEElectronCrossSectionModel(const G4String& sss)
{
  if(IsLocked()) { return; }
  fCParameters->SetPIXEElectronCrossSectionModel(sss);
}

const G4String& G4EmParameters::PIXEElectronCrossSectionModel()
{
  return fCParameters->PIXEElectronCrossSectionModel();
}

void G4EmParameters::SetLivermoreDataDir(const G4String& sss)
{
  if(IsLocked()) { return; }
  fCParameters->SetLivermoreDataDir(sss);
}

const G4String& G4EmParameters::LivermoreDataDir()
{
  return fCParameters->LivermoreDataDir();
}

void G4EmParameters::PrintWarning(G4ExceptionDescription& ed) const
{
  G4Exception("G4EmParameters", "em0044", JustWarning, ed);
}

void G4EmParameters::AddPAIModel(const G4String& particle,
                                 const G4String& region,
                                 const G4String& type)
{
  if(IsLocked()) { return; }
  fBParameters->AddPAIModel(particle, region, type);
}

const std::vector<G4String>& G4EmParameters::ParticlesPAI() const
{
  return fBParameters->ParticlesPAI();
}

const std::vector<G4String>& G4EmParameters::RegionsPAI() const
{
  return fBParameters->RegionsPAI();
}

const std::vector<G4String>& G4EmParameters::TypesPAI() const
{
  return fBParameters->TypesPAI();
}

void G4EmParameters::AddMicroElec(const G4String& region)
{
  if(IsLocked()) { return; }
  fCParameters->AddMicroElec(region);
}

const std::vector<G4String>& G4EmParameters::RegionsMicroElec() const
{
  return fCParameters->RegionsMicroElec();
}

void G4EmParameters::AddDNA(const G4String& region, const G4String& type)
{
  if(IsLocked()) { return; }
  fCParameters->AddDNA(region, type);
  ActivateDNA();
}

const std::vector<G4String>& G4EmParameters::RegionsDNA() const
{
  return fCParameters->RegionsDNA();
}

const std::vector<G4String>& G4EmParameters::TypesDNA() const
{
  return fCParameters->TypesDNA();
}

void G4EmParameters::AddPhysics(const G4String& region, const G4String& type)
{
  if(IsLocked()) { return; }
  fBParameters->AddPhysics(region, type);
}

const std::vector<G4String>& G4EmParameters::RegionsPhysics() const
{
  return fBParameters->RegionsPhysics();
}

const std::vector<G4String>& G4EmParameters::TypesPhysics() const
{
  return fBParameters->TypesPhysics();
}

void G4EmParameters::SetSubCutRegion(const G4String& region)
{
  if(IsLocked()) { return; }
  fBParameters->SetSubCutRegion(region);
}

void 
G4EmParameters::SetDeexActiveRegion(const G4String& region, G4bool adeex,
                                    G4bool aauger, G4bool apixe)
{
  if(IsLocked()) { return; }
  fCParameters->SetDeexActiveRegion(region, adeex, aauger, apixe);
}

void 
G4EmParameters::SetProcessBiasingFactor(const G4String& procname, 
                                        G4double val, G4bool wflag)
{
  if(IsLocked()) { return; }
  fBParameters->SetProcessBiasingFactor(procname, val, wflag);
}

void 
G4EmParameters::ActivateForcedInteraction(const G4String& procname, 
                                          const G4String& region,
                                          G4double length, 
                                          G4bool wflag)
{
  if(IsLocked() && !gener) { return; }
  fBParameters->ActivateForcedInteraction(procname, region, length, wflag);
}

void 
G4EmParameters::ActivateSecondaryBiasing(const G4String& procname,
                                         const G4String& region, 
                                         G4double factor,
                                         G4double energyLim)
{
  if(IsLocked()) { return; }
  fBParameters->ActivateSecondaryBiasing(procname, region, factor, energyLim);
}

void G4EmParameters::DefineRegParamForLoss(G4VEnergyLossProcess* ptr) const
{
  fBParameters->DefineRegParamForLoss(ptr);
}

void G4EmParameters::DefineRegParamForEM(G4VEmProcess* ptr) const
{
  fBParameters->DefineRegParamForEM(ptr);
}

G4bool G4EmParameters::QuantumEntanglement() const
{
  return fBParameters->QuantumEntanglement(); 
}

void G4EmParameters::SetQuantumEntanglement(G4bool v)
{
  if(IsLocked()) { return; }
  fBParameters->SetQuantumEntanglement(v); 
}

G4bool G4EmParameters::GetDirectionalSplitting() const { 
  return fBParameters->GetDirectionalSplitting(); 
}

void G4EmParameters::SetDirectionalSplitting(G4bool v) 
{ 
  if(IsLocked()) { return; }
  fBParameters->SetDirectionalSplitting(v); 
}

void G4EmParameters::SetDirectionalSplittingTarget(const G4ThreeVector& v)
{ 
  if(IsLocked()) { return; }
  fBParameters->SetDirectionalSplittingTarget(v);
}

G4ThreeVector G4EmParameters::GetDirectionalSplittingTarget() const
{ 
  return fBParameters->GetDirectionalSplittingTarget(); 
}

void G4EmParameters::SetDirectionalSplittingRadius(G4double r)
{ 
  if(IsLocked()) { return; }
  fBParameters->SetDirectionalSplittingRadius(r); 
}

G4double G4EmParameters::GetDirectionalSplittingRadius()
{ 
  return fBParameters->GetDirectionalSplittingRadius(); 
}

void G4EmParameters::DefineRegParamForDeex(G4VAtomDeexcitation* ptr) const
{
  fCParameters->DefineRegParamForDeex(ptr); 
}

void G4EmParameters::StreamInfo(std::ostream& os) const
{
  G4long prec = os.precision(5);
  os << "=======================================================================" << "\n";
  os << "======                 Electromagnetic Physics Parameters      ========" << "\n";
  os << "=======================================================================" << "\n";
  os << "LPM effect enabled                                 " <<flagLPM << "\n";
  os << "Enable creation and use of sampling tables         " <<fSamplingTable << "\n";
  os << "Apply cuts on all EM processes                     " <<applyCuts << "\n";
  const char* transportationWithMsc = "Disabled";
  if(fTransportationWithMsc == G4TransportationWithMscType::fEnabled) {
    transportationWithMsc = "Enabled";
  } else if (fTransportationWithMsc == G4TransportationWithMscType::fMultipleSteps) {
    transportationWithMsc = "MultipleSteps";
  }
  os << "Use combined TransportationWithMsc                 " <<transportationWithMsc << "\n";
  os << "Use general process                                " <<gener << "\n";
  os << "Enable linear polarisation for gamma               " <<fPolarisation << "\n";
  os << "Enable photoeffect sampling below K-shell          " <<fPEKShell << "\n";
  os << "Enable sampling of quantum entanglement            " 
     <<fBParameters->QuantumEntanglement()  << "\n";
  os << "X-section factor for integral approach             " <<lambdaFactor << "\n";
  os << "Min kinetic energy for tables                      " 
     <<G4BestUnit(minKinEnergy,"Energy") << "\n";
  os << "Max kinetic energy for tables                      " 
     <<G4BestUnit(maxKinEnergy,"Energy") << "\n";
  os << "Number of bins per decade of a table               " <<nbinsPerDecade << "\n";
  os << "Verbose level                                      " <<verbose << "\n";
  os << "Verbose level for worker thread                    " <<workerVerbose << "\n";
  os << "Bremsstrahlung energy threshold above which \n" 
     << "  primary e+- is added to the list of secondary    " 
     <<G4BestUnit(bremsTh,"Energy") << "\n";
  os << "Bremsstrahlung energy threshold above which primary\n" 
     << "  muon/hadron is added to the list of secondary    " 
     <<G4BestUnit(bremsMuHadTh,"Energy") << "\n";
  os << "Lowest triplet kinetic energy                      " 
     <<G4BestUnit(lowestTripletEnergy,"Energy") << "\n";
  os << "Enable sampling of gamma linear polarisation       " <<fPolarisation << "\n";
  os << "5D gamma conversion model type                     " <<tripletConv << "\n";
  os << "5D gamma conversion model on isolated ion          " <<onIsolated << "\n";
  if(max5DEnergyForMuPair>0.0) {
  os << "5D gamma conversion limit for muon pair            " 
     << max5DEnergyForMuPair/CLHEP::GeV << " GeV\n";
  }
  os << "Livermore data directory                           " 
     << fCParameters->LivermoreDataDir() << "\n";

  os << "=======================================================================" << "\n";
  os << "======                 Ionisation Parameters                   ========" << "\n";
  os << "=======================================================================" << "\n";
  os << "Step function for e+-                              " 
     <<"("<<fBParameters->GetStepFunctionP1() << ", " 
     << fBParameters->GetStepFunctionP2()/CLHEP::mm << " mm)\n";
  os << "Step function for muons/hadrons                    " 
     <<"("<<fBParameters->GetStepFunctionMuHadP1() << ", " 
     << fBParameters->GetStepFunctionMuHadP2()/CLHEP::mm << " mm)\n";
  os << "Step function for light ions                       " 
     <<"("<<fBParameters->GetStepFunctionLightIonsP1() << ", " 
     << fBParameters->GetStepFunctionLightIonsP2()/CLHEP::mm << " mm)\n";
  os << "Step function for general ions                     " 
     <<"("<<fBParameters->GetStepFunctionIonsP1() << ", " 
     << fBParameters->GetStepFunctionIonsP2()/CLHEP::mm << " mm)\n";
  os << "Lowest e+e- kinetic energy                         " 
     <<G4BestUnit(lowestElectronEnergy,"Energy") << "\n";
  os << "Lowest muon/hadron kinetic energy                  " 
     <<G4BestUnit(lowestMuHadEnergy,"Energy") << "\n";
  os << "Use ICRU90 data                                    " << fICRU90 << "\n";
  os << "Fluctuations of dE/dx are enabled                  " <<lossFluctuation << "\n";
  G4String namef = "Universal";
  if(fFluct == fUrbanFluctuation) { namef = "Urban"; }
  else if(fFluct == fDummyFluctuation) { namef = "Dummy"; }
  os << "Type of fluctuation model for leptons and hadrons  " << namef << "\n";
  os << "Use built-in Birks satuaration                     " << birks << "\n";
  os << "Build CSDA range enabled                           " <<buildCSDARange << "\n";
  os << "Use cut as a final range enabled                   " <<cutAsFinalRange << "\n";
  os << "Enable angular generator interface                 " 
     <<useAngGeneratorForIonisation << "\n";
  os << "Max kinetic energy for CSDA tables                 " 
     <<G4BestUnit(maxKinEnergyCSDA,"Energy") << "\n";
  os << "Max kinetic energy for NIEL computation            " 
     <<G4BestUnit(maxNIELEnergy,"Energy") << "\n";
  os << "Linear loss limit                                  " <<linLossLimit << "\n";
  os << "Read data from file for e+e- pair production by mu " <<fMuDataFromFile << "\n";

  os << "=======================================================================" << "\n";
  os << "======                 Multiple Scattering Parameters          ========" << "\n";
  os << "=======================================================================" << "\n";
  os << "Type of msc step limit algorithm for e+-           " <<mscStepLimit << "\n";
  os << "Type of msc step limit algorithm for muons/hadrons " <<mscStepLimitMuHad << "\n";
  os << "Msc lateral displacement for e+- enabled           " <<lateralDisplacement << "\n";
  os << "Msc lateral displacement for muons and hadrons     " <<muhadLateralDisplacement << "\n";
  os << "Urban msc model lateral displacement alg96         " <<lateralDisplacementAlg96 << "\n";
  os << "Range factor for msc step limit for e+-            " <<rangeFactor << "\n";
  os << "Range factor for msc step limit for muons/hadrons  " <<rangeFactorMuHad << "\n";
  os << "Geometry factor for msc step limitation of e+-     " <<geomFactor << "\n";
  os << "Safety factor for msc step limit for e+-           " <<safetyFactor << "\n";
  os << "Skin parameter for msc step limitation of e+-      " <<skin << "\n";
  os << "Lambda limit for msc step limit for e+-            " <<lambdaLimit/CLHEP::mm << " mm\n";
  os << "Use Mott correction for e- scattering              " << useMottCorrection << "\n";
  os << "Factor used for dynamic computation of angular \n" 
     << "  limit between single and multiple scattering     " << factorForAngleLimit << "\n";
  os << "Fixed angular limit between single \n"
     << "  and multiple scattering                          " 
     << thetaLimit/CLHEP::rad << " rad\n";
  os << "Upper energy limit for e+- multiple scattering     " 
     << energyLimit/CLHEP::MeV << " MeV\n";
  os << "Type of electron single scattering model           " <<fSStype << "\n";
  os << "Type of nuclear form-factor                        " <<nucFormfactor << "\n";
  os << "Screening factor                                   " <<factorScreen << "\n";
  os << "=======================================================================" << "\n";

  if(fCParameters->Fluo()) {
  os << "======                 Atomic Deexcitation Parameters          ========" << "\n";
  os << "=======================================================================" << "\n";
  os << "Fluorescence enabled                               " <<fCParameters->Fluo() << "\n";
  G4String named = "fluor";
  G4EmFluoDirectory fdir = FluoDirectory();
  if(fdir == fluoBearden) { named = "fluor_Bearden"; }
  else if(fdir == fluoANSTO) { named = "fluor_ANSTO"; }
  else if(fdir == fluoXDB_EADL) { named = "fluor_XDB_EADL"; }
  os << "Directory in G4LEDATA for fluorescence data files  " << named << "\n";
  os << "Auger electron cascade enabled                     " 
     <<fCParameters->Auger() << "\n";
  os << "PIXE atomic de-excitation enabled                  " <<fCParameters->Pixe() << "\n";
  os << "De-excitation module ignores cuts                  " 
     <<fCParameters->DeexcitationIgnoreCut() << "\n";
  os << "Type of PIXE cross section for hadrons             " 
     <<fCParameters->PIXECrossSectionModel() << "\n";
  os << "Type of PIXE cross section for e+-                 " 
     <<fCParameters->PIXEElectronCrossSectionModel() << "\n";
  os << "=======================================================================" << "\n";
  }
  if(fDNA) {
  os << "======                 DNA Physics Parameters                  ========" << "\n";
  os << "=======================================================================" << "\n";
  os << "Use fast sampling in DNA models                    " 
     << fCParameters->DNAFast() << "\n";
  os << "Use Stationary option in DNA models                " 
     << fCParameters->DNAStationary() << "\n";
  os << "Use DNA with multiple scattering of e-             " 
     << fCParameters->DNAElectronMsc() << "\n";
  os << "Use DNA e- solvation model type                    " 
     << fCParameters->DNAeSolvationSubType() << "\n";
  os << "=======================================================================" << G4endl;
  }
  os.precision(prec);
}

void G4EmParameters::Dump()
{
  if(fIsPrinted) return;

#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&emParametersMutex);
#endif
  StreamInfo(G4cout);
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&emParametersMutex);
#endif
}

std::ostream& operator<< (std::ostream& os, const G4EmParameters& par)
{
  par.StreamInfo(os);
  return os;
}

G4bool G4EmParameters::IsLocked() const
{
  return (!G4Threading::IsMasterThread() ||
	  (fStateManager->GetCurrentState() != G4State_PreInit &&
       	   fStateManager->GetCurrentState() != G4State_Init &&
	   fStateManager->GetCurrentState() != G4State_Idle));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
