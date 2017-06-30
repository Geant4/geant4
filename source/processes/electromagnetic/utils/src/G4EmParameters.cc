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
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VEmProcess.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4EmParametersMessenger.hh"
#include "G4NistManager.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ApplicationState.hh"
#include "G4StateManager.hh"
#include "G4Threading.hh"

G4EmParameters* G4EmParameters::theInstance = nullptr;

#ifdef G4MULTITHREADED
  G4Mutex G4EmParameters::emParametersMutex = G4MUTEX_INITIALIZER;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmParameters* G4EmParameters::Instance()
{
  if(nullptr == theInstance) {
    static G4EmParameters manager;
    theInstance = &manager;
  }
  return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmParameters::~G4EmParameters()
{
  delete theMessenger;
  delete emSaturation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmParameters::G4EmParameters()
{
  G4NistManager::Instance();
  theMessenger = new G4EmParametersMessenger(this);

  fStateManager = G4StateManager::GetStateManager();
  Initialise();
  emSaturation = nullptr;
}

void G4EmParameters::SetDefaults()
{
  if(!IsLocked()) { Initialise(); }
}

void G4EmParameters::Initialise()
{
  lossFluctuation = true;
  buildCSDARange = false;
  flagLPM = true;
  spline = true;
  cutAsFinalRange = false;
  applyCuts = false;
  fluo = false;
  beardenFluoDir = false;
  auger = false;
  augerCascade = false;
  pixe = false;
  deexIgnoreCut = false;
  lateralDisplacement = true;
  muhadLateralDisplacement = false;
  latDisplacementBeyondSafety = false;
  useAngGeneratorForIonisation = false;
  useMottCorrection = false;
  integral = true;
  birks = false;

  minSubRange = 1.0;
  minKinEnergy = 0.1*CLHEP::keV;
  maxKinEnergy = 100.0*CLHEP::TeV;
  maxKinEnergyCSDA = 1.0*CLHEP::GeV;
  lowestElectronEnergy = 1.0*CLHEP::keV;
  lowestMuHadEnergy = 1.0*CLHEP::keV;
  linLossLimit = 0.01;
  bremsTh = maxKinEnergy;
  lambdaFactor = 0.8;
  factorForAngleLimit = 1.0;
  thetaLimit = CLHEP::pi;
  rangeFactor = 0.04;
  rangeFactorMuHad = 0.2;
  geomFactor = 2.5;
  skin = 1.0;
  dRoverRange = 0.2;
  finalRange = CLHEP::mm;
  dRoverRangeMuHad = 0.2;
  finalRangeMuHad = 0.1*CLHEP::mm;
  factorScreen = 1.0;

  nbins  = 77;
  nbinsPerDecade = 7;
  verbose = 1;
  workerVerbose = 0;

  mscStepLimit = fUseSafety;
  mscStepLimitMuHad = fMinimal;
  nucFormfactor = fExponentialNF;

  namePIXE = "Empirical";
  nameElectronPIXE = "Livermore";
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

void G4EmParameters::SetSpline(G4bool val)
{
  if(IsLocked()) { return; }
  spline = val;
}

G4bool G4EmParameters::Spline() const
{
  return spline;
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
  fluo = val;
}

G4bool G4EmParameters::Fluo() const
{
  return fluo;
}

void G4EmParameters::SetBeardenFluoDir(G4bool val)
{
  if(IsLocked()) { return; }
  beardenFluoDir = val;
}

G4bool G4EmParameters::BeardenFluoDir() const
{
  return beardenFluoDir;
}

void G4EmParameters::SetAuger(G4bool val)
{
  if(IsLocked()) { return; }
  auger = val;
  if(val) { fluo = true; }
}

G4bool G4EmParameters::Auger() const
{
  return auger;
}

void G4EmParameters::SetAugerCascade(G4bool val)
{
  if(IsLocked()) { return; }
  augerCascade = val;
  if(val) { fluo = true; auger = true; }
}

G4bool G4EmParameters::AugerCascade() const
{
  return augerCascade;
}

void G4EmParameters::SetPixe(G4bool val)
{
  if(IsLocked()) { return; }
  pixe = val;
  if(val) { fluo = true; }
}

G4bool G4EmParameters::Pixe() const
{
  return pixe;
}

void G4EmParameters::SetDeexcitationIgnoreCut(G4bool val)
{
  if(IsLocked()) { return; }
  deexIgnoreCut = val;
}

G4bool G4EmParameters::DeexcitationIgnoreCut() const
{
  return deexIgnoreCut;
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

void G4EmParameters::SetMuHadLateralDisplacement(G4bool val)
{
  if(IsLocked()) { return; }
  muhadLateralDisplacement = val;
}

G4bool G4EmParameters::MuHadLateralDisplacement() const
{
  return muhadLateralDisplacement;
}

void G4EmParameters::SetLatDisplacementBeyondSafety(G4bool val)
{
  if(IsLocked()) { return; }
  latDisplacementBeyondSafety = val;
}

G4bool G4EmParameters::LatDisplacementBeyondSafety() const
{
  return latDisplacementBeyondSafety;
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

void G4EmParameters::SetBirksActive(G4bool val)
{
  birks = val;
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&G4EmParameters::emParametersMutex);
#endif
  if(birks) {
    if(!emSaturation) { emSaturation = new G4EmSaturation(1); }
    emSaturation->InitialiseG4Saturation();
  }
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&G4EmParameters::emParametersMutex);
#endif
}

G4bool G4EmParameters::BirksActive() const
{
  return birks;
}

void G4EmParameters::SetEmSaturation(G4EmSaturation* ptr)
{
  if(emSaturation != ptr) {
    delete emSaturation;
    emSaturation = ptr;
    SetBirksActive(true);
  }
}

G4EmSaturation* G4EmParameters::GetEmSaturation()
{
  if(!emSaturation) { SetBirksActive(true); }
  return emSaturation;
}

void G4EmParameters::SetMinSubRange(G4double val)
{
  if(IsLocked()) { return; }
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
  if(IsLocked()) { return; }
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
  if(IsLocked()) { return; }
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
  if(IsLocked()) { return; }
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

void G4EmParameters::SetLowestElectronEnergy(G4double val)
{
  if(IsLocked()) { return; }
  if(val >= 0.0) {
    lowestElectronEnergy = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of lowestElectronEnergy is out of range: " 
       << val/MeV << " MeV is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::LowestElectronEnergy() const
{
  return lowestElectronEnergy; 
}

void G4EmParameters::SetLowestMuHadEnergy(G4double val)
{
  if(IsLocked()) { return; }
  if(val >= 0.0) {
    lowestMuHadEnergy = val;
  } else {
    G4ExceptionDescription ed;
    ed << "Value of lowestMuHadEnergy is out of range: " 
       << val/MeV << " MeV is ignored"; 
    PrintWarning(ed);
  }
}

G4double G4EmParameters::LowestMuHadEnergy() const
{
  return lowestMuHadEnergy; 
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
  if(v1 > 0.0 && v1 <= 1.0 && v2 > 0.0) {
    dRoverRange = v1;
    finalRange = v2;
  } else {
    G4ExceptionDescription ed;
    ed << "Values of step function are out of range: " 
       << v1 << ", " << v2/CLHEP::mm << " mm - are ignored"; 
    PrintWarning(ed);
  }
}

void G4EmParameters::SetStepFunctionMuHad(G4double v1, G4double v2)
{
  if(IsLocked()) { return; }
  if(v1 > 0.0 && v1 <= 1.0 && v2 > 0.0) {
    dRoverRangeMuHad = v1;
    finalRangeMuHad = v2;
  } else {
    G4ExceptionDescription ed;
    ed << "Values of step function are out of range: " 
       << v1 << ", " << v2/CLHEP::mm << " mm - are ignored"; 
    PrintWarning(ed);
  }
}

void G4EmParameters::SetNumberOfBins(G4int val)
{
  if(IsLocked()) { return; }
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
  if(IsLocked()) { return; }
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

void G4EmParameters::SetPIXECrossSectionModel(const G4String& sss)
{
  G4cout << "G4EmParameters::SetPIXECrossSectionModel " << sss << G4endl;
  if(IsLocked()) { return; }
  namePIXE = sss;
}

const G4String& G4EmParameters::PIXECrossSectionModel()
{
  return namePIXE;
}

void G4EmParameters::SetPIXEElectronCrossSectionModel(const G4String& sss)
{
  if(IsLocked()) { return; }
  nameElectronPIXE = sss;
}

const G4String& G4EmParameters::PIXEElectronCrossSectionModel()
{
  return nameElectronPIXE;
}

void G4EmParameters::PrintWarning(G4ExceptionDescription& ed) const
{
  G4Exception("G4EmParameters", "em0044", JustWarning, ed);
}

G4String G4EmParameters::CheckRegion(const G4String& reg) const
{
  G4String r = reg;
  if(r == "" || r == "world" || r == "World") {
    r = "DefaultRegionForTheWorld";
  }
  return r;
}

void G4EmParameters::AddPAIModel(const G4String& particle,
                                 const G4String& region,
                                 const G4String& type)
{
  G4String r = CheckRegion(region);
  G4int nreg =  m_regnamesPAI.size();
  for(G4int i=0; i<nreg; ++i) {
    if((m_particlesPAI[i] == particle || 
        m_particlesPAI[i] == "all" || 
        particle == "all") && 
       (m_regnamesPAI[i] == r || 
        m_regnamesPAI[i] == "DefaultRegionForTheWorld" || 
        r == "DefaultRegionForTheWorld") ) {

      m_typesPAI[i] = type;
      if(particle == "all") { m_particlesPAI[i] = particle; }
      if(r == "DefaultRegionForTheWorld") { m_regnamesPAI[i] = r; }
      return;
    }
  }
  m_particlesPAI.push_back(particle);
  m_regnamesPAI.push_back(r);
  m_typesPAI.push_back(type);
}

const std::vector<G4String>& G4EmParameters::ParticlesPAI() const
{
  return m_particlesPAI;
}

const std::vector<G4String>& G4EmParameters::RegionsPAI() const
{
  return m_regnamesPAI;
}

const std::vector<G4String>& G4EmParameters::TypesPAI() const
{
  return m_typesPAI;
}

void G4EmParameters::AddMicroElec(const G4String& region)
{
  G4String r = CheckRegion(region);
  G4int nreg =  m_regnamesME.size();
  for(G4int i=0; i<nreg; ++i) {
    if(r == m_regnamesME[i]) { return; }
  }
  m_regnamesME.push_back(r);
}

const std::vector<G4String>& G4EmParameters::RegionsMicroElec() const
{
  return m_regnamesME;
}

void G4EmParameters::AddDNA(const G4String& region, const G4String& type)
{
  G4String r = CheckRegion(region);
  G4int nreg =  m_regnamesDNA.size();
  for(G4int i=0; i<nreg; ++i) {
    if(r == m_regnamesDNA[i]) { return; }
  }
  m_regnamesDNA.push_back(r);
  m_typesDNA.push_back(type);
}

const std::vector<G4String>& G4EmParameters::RegionsDNA() const
{
  return m_regnamesDNA;
}

const std::vector<G4String>& G4EmParameters::TypesDNA() const
{
  return m_typesDNA;
}

void G4EmParameters::AddMsc(const G4String& region, const G4String& type)
{
  G4String r = CheckRegion(region);
  G4int nreg =  m_regnamesMsc.size();
  for(G4int i=0; i<nreg; ++i) {
    if(r == m_regnamesMsc[i]) { return; }
  }
  m_regnamesMsc.push_back(r);
  m_typesMsc.push_back(type);
}

const std::vector<G4String>& G4EmParameters::RegionsMsc() const
{
  return m_regnamesMsc;
}

const std::vector<G4String>& G4EmParameters::TypesMsc() const
{
  return m_typesMsc;
}

void G4EmParameters::AddPhysics(const G4String& region, const G4String& type)
{
  G4String r = CheckRegion(region);
  G4int nreg =  m_regnamesMsc.size();
  for(G4int i=0; i<nreg; ++i) {
    if(r == m_regnamesMsc[i]) { return; }
  }
  m_regnamesMsc.push_back(r);
  m_typesMsc.push_back(type);
}

const std::vector<G4String>& G4EmParameters::RegionsPhysics() const
{
  return m_regnamesMsc;
}

const std::vector<G4String>& G4EmParameters::TypesPhysics() const
{
  return m_typesMsc;
}

void G4EmParameters::SetSubCutoff(G4bool val, const G4String& region)
{
  if(IsLocked()) { return; }
  G4String r = CheckRegion(region);
  G4int nreg =  m_regnamesSubCut.size();
  for(G4int i=0; i<nreg; ++i) {
    if(r == m_regnamesSubCut[i]) { 
      m_subCuts[i] = val;
      return; 
    }
  }
  m_regnamesSubCut.push_back(r);
  m_subCuts.push_back(val);
}

void 
G4EmParameters::SetDeexActiveRegion(const G4String& region, G4bool fdeex,
                                    G4bool fauger, G4bool fpixe)
{
  if(IsLocked()) { return; }
  if(fdeex) { fluo = true; }
  G4String r = CheckRegion(region);
  G4int nreg =  m_regnamesDeex.size();
  if(0 == nreg && r != "DefaultRegionForTheWorld") {
    m_regnamesDeex.push_back("DefaultRegionForTheWorld");
    m_fluo.push_back(false);
    m_auger.push_back(false);
    m_pixe.push_back(false);
    nreg = 1;
  }
  for(G4int i=0; i<nreg; ++i) {
    if(r == m_regnamesDeex[i]) { 
      m_fluo[i] = fdeex;
      m_auger[i]= fauger;
      m_pixe[i] = fpixe;
      return; 
    }
  }
  m_regnamesDeex.push_back(r);
  m_fluo.push_back(fdeex);
  m_auger.push_back(fauger);
  m_pixe.push_back(fpixe);
}

void 
G4EmParameters::SetProcessBiasingFactor(const G4String& procname, 
                                        G4double val, G4bool wflag)
{
  if(IsLocked()) { return; }
  if(val > 0.0) {
    G4int n =  m_procBiasedXS.size();
    for(G4int i=0; i<n; ++i) {
      if(procname == m_procBiasedXS[i]) { 
	m_factBiasedXS[i] = val;
	m_weightBiasedXS[i]= wflag;
	return; 
      }
    }
    m_procBiasedXS.push_back(procname);
    m_factBiasedXS.push_back(val);
    m_weightBiasedXS.push_back(wflag);
  } else {
    G4ExceptionDescription ed;
    ed << "Process: " << procname << " XS biasing factor " 
       << val << " is negative - ignored"; 
    PrintWarning(ed);
  }
}

void 
G4EmParameters::ActivateForcedInteraction(const G4String& procname, 
                                          const G4String& region,
                                          G4double length, 
                                          G4bool wflag)
{
  if(IsLocked()) { return; }
  G4String r = CheckRegion(region);
  if(length >= 0.0) {
    G4int n =  m_procForced.size();
    for(G4int i=0; i<n; ++i) {
      if(procname == m_procForced[i] && r == m_regnamesForced[i] ) { 
	m_lengthForced[i] = length;
	m_weightForced[i]= wflag;
	return; 
      }
    }
    m_regnamesForced.push_back(r);
    m_procForced.push_back(procname);
    m_lengthForced.push_back(length);
    m_weightForced.push_back(wflag);
  } else {
    G4ExceptionDescription ed;
    ed << "Process: " << procname << " in region " << r
       << " : forced interacttion length= " 
       << length << " is negative - ignored"; 
    PrintWarning(ed);
  }
}

void 
G4EmParameters::ActivateSecondaryBiasing(const G4String& procname,
                                         const G4String& region, 
                                         G4double factor,
                                         G4double energyLimit)
{
  if(IsLocked()) { return; }
  G4String r = CheckRegion(region);
  if(factor >= 0.0 && energyLimit >= 0.0) {
    G4int n =  m_procBiasedSec.size();
    for(G4int i=0; i<n; ++i) {
      if(procname == m_procBiasedSec[i] && r == m_regnamesBiasedSec[i] ) { 
	m_factBiasedSec[i] = factor;
	m_elimBiasedSec[i] = energyLimit;
	return; 
      }
    }
    m_regnamesBiasedSec.push_back(r);
    m_procBiasedSec.push_back(procname);
    m_factBiasedSec.push_back(factor);
    m_elimBiasedSec.push_back(energyLimit);
  } else {
    G4ExceptionDescription ed;
    ed << "Process: " << procname << " in region " << r
       << " : secondary bised factor= " 
       << factor << ", Elim= " << energyLimit <<  " - ignored"; 
    PrintWarning(ed);
  }
}

void G4EmParameters::DefineRegParamForLoss(G4VEnergyLossProcess* ptr, 
                                           G4bool isElectron) const
{
  if(isElectron) { ptr->SetStepFunction(dRoverRange, finalRange, false); }
  else { ptr->SetStepFunction(dRoverRangeMuHad, finalRangeMuHad, false); }

  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  G4int n = m_regnamesSubCut.size();
  for(G4int i=0; i<n; ++i) { 
    const G4Region* reg = regionStore->GetRegion(m_regnamesSubCut[i], false);
    if(reg) { ptr->ActivateSubCutoff(m_subCuts[i], reg); }
  }
  n = m_procBiasedXS.size();
  for(G4int i=0; i<n; ++i) {
    if(ptr->GetProcessName() == m_procBiasedXS[i]) {
      ptr->SetCrossSectionBiasingFactor(m_factBiasedXS[i], 
					m_weightBiasedXS[i]);
      break; 
    }
  }
  n = m_procForced.size();
  for(G4int i=0; i<n; ++i) {
    if(ptr->GetProcessName() == m_procForced[i]) {
      ptr->ActivateForcedInteraction(m_lengthForced[i],
				     m_regnamesForced[i],
				     m_weightForced[i]);
      break; 
    }
  }
  n = m_procBiasedSec.size();
  for(G4int i=0; i<n; ++i) {
    if(ptr->GetProcessName() == m_procBiasedSec[i]) {
      ptr->ActivateSecondaryBiasing(m_regnamesBiasedSec[i],
				    m_factBiasedSec[i], 
				    m_elimBiasedSec[i]);
      break; 
    }
  }
}

void G4EmParameters::DefineRegParamForEM(G4VEmProcess* ptr) const
{
  G4int n = m_procBiasedXS.size();
  for(G4int i=0; i<n; ++i) {
    if(ptr->GetProcessName() == m_procBiasedXS[i]) {
      ptr->SetCrossSectionBiasingFactor(m_factBiasedXS[i], 
					m_weightBiasedXS[i]);
      break; 
    }
  }
  n = m_procForced.size();
  for(G4int i=0; i<n; ++i) {
    if(ptr->GetProcessName() == m_procForced[i]) {
      ptr->ActivateForcedInteraction(m_lengthForced[i],
				     m_regnamesForced[i],
				     m_weightForced[i]);
      break; 
    }
  }
  n = m_procBiasedSec.size();
  for(G4int i=0; i<n; ++i) {
    if(ptr->GetProcessName() == m_procBiasedSec[i]) {
      ptr->ActivateSecondaryBiasing(m_regnamesBiasedSec[i],
				    m_factBiasedSec[i], 
				    m_elimBiasedSec[i]);
      break; 
    }
  }
}

void G4EmParameters::DefineRegParamForDeex(G4VAtomDeexcitation* ptr) const
{
  G4int n = m_regnamesDeex.size();
  for(G4int i=0; i<n; ++i) {
    ptr->SetDeexcitationActiveRegion(m_regnamesDeex[i],
				     m_fluo[i], m_auger[i], m_pixe[i]);
  }
}

std::ostream& G4EmParameters::StreamInfo(std::ostream& os) const
{
  G4int prec = os.precision(5);
  os << "=======================================================================" << "\n";
  os << "======                 Electromagnetic Physics Parameters      ========" << "\n";
  os << "=======================================================================" << "\n";
  os << "Fluctuations of dE/dx are enabled                  " <<lossFluctuation << "\n";
  os << "Build CSDA range enabled                           " <<buildCSDARange << "\n";
  os << "LPM effect enabled                                 " <<flagLPM << "\n";
  os << "Spline of EM tables enabled                        " <<spline << "\n";
  os << "Use cut as a final range enabled                   " <<finalRange << "\n";
  os << "Apply cuts on all EM processes                     " <<applyCuts << "\n";
  os << "Fluorescence enabled                               " <<fluo << "\n";
  os << "Fluorescence Bearden data files enabled            " <<beardenFluoDir << "\n";
  os << "Auger electron production enabled                  " <<auger << "\n";
  os << "Auger cascade enabled                              " <<augerCascade << "\n";
  os << "PIXE atomic de-excitation enabled                  " <<pixe << "\n";
  os << "De-excitation module ignores cuts                  " <<deexIgnoreCut << "\n";
  os << "Msc lateraral displacement for e+- enabled         " <<lateralDisplacement << "\n";
  os << "Msc lateraral displacement for muons and hadrons   " <<muhadLateralDisplacement << "\n";
  os << "Msc lateraral displacement beyond geometry safety  " <<latDisplacementBeyondSafety << "\n";
  os << "Enable angular generator interface                 " 
     <<useAngGeneratorForIonisation << "\n";
  os << "Use Mott correction for e- scattering              " 
     <<useMottCorrection << "\n";
  os << "Use integral approach for tracking                 " 
     <<integral << "\n";
  os << "Use built-in Birks satuaration                     " 
     << birks << "\n";

  os << "Factor of cut reduction for sub-cutoff method      " <<minSubRange << "\n";
  os << "Min kinetic energy for tables                      " 
     <<G4BestUnit(minKinEnergy,"Energy") << "\n";
  os << "Max kinetic energy for tables                      " 
     <<G4BestUnit(maxKinEnergy,"Energy") << "\n";
  os << "Max kinetic energy for CSDA tables                 " 
     <<G4BestUnit(maxKinEnergyCSDA,"Energy") << "\n";
  os << "Lowest e+e- kinetic energy                         " 
     <<G4BestUnit(lowestElectronEnergy,"Energy") << "\n";
  os << "Lowest muon/hadron kinetic energy                  " 
     <<G4BestUnit(lowestMuHadEnergy,"Energy") << "\n";
  os << "Linear loss limit " <<linLossLimit << "\n";
  os << "Bremsstrahlung energy threshold above which \n" 
     << "  primary is added to the list of secondary        " 
     <<G4BestUnit(bremsTh,"Energy") << "\n";
  os << "X-section factor for integral approach             " <<lambdaFactor << "\n";
  os << "Factor used for dynamic computation of angular \n" 
     << "  limit between single and multiple scattering     " << factorForAngleLimit << "\n";
  os << "Fixed angular limit between single \n"
     << "  and multiple scattering                          " 
     <<thetaLimit/rad << " rad" << "\n";
  os << "Range factor for msc step limit for e+-            " <<rangeFactor << "\n";
  os << "Range factor for msc step limit for muons/hadrons  " <<rangeFactorMuHad << "\n";
  os << "Geometry factor for msc step limitation of e+-     " <<geomFactor << "\n";
  os << "Skin parameter for msc step limitation of e+-      " <<skin << "\n";
  os << "Screening factor                                   " <<factorScreen << "\n";
  os << "Step function for e+-                              " <<"("<< dRoverRange
     << ", " << finalRange << " mm)\n";
  os << "Step function for muons/hadrons                    " <<"("<< dRoverRangeMuHad
     << ", " << finalRangeMuHad << " mm)\n";

  os << "Number of bins in tables                           " <<nbins   << "\n";
  os << "Number of bins per decade of a table               " <<nbinsPerDecade << "\n";
  os << "Verbose level                                      " <<verbose << "\n";
  os << "Verbose level for worker thread                    " <<workerVerbose << "\n";

  os << "Type of msc step limit algorithm for e+-           " <<mscStepLimit << "\n";
  os << "Type of msc step limit algorithm for muons/hadrons " <<mscStepLimitMuHad << "\n";
  os << "Type of nuclear form-factor                        " <<nucFormfactor << "\n";

  os << "Type of PIXE cross section for hadrons             " <<namePIXE << "\n";
  os << "Type of PIXE cross section for e+-                 " <<nameElectronPIXE << "\n";
  os << "=======================================================================" << "\n";
  os.precision(prec);
  return os;
}

void G4EmParameters::Dump() const
{
  StreamInfo(G4cout);
}

std::ostream& operator<< (std::ostream& os, const G4EmParameters& par)
{
  return par.StreamInfo(os);
}

G4bool G4EmParameters::IsLocked() const
{
  return (!G4Threading::IsMasterThread() ||
	  (fStateManager->GetCurrentState() != G4State_PreInit &&
	   fStateManager->GetCurrentState() != G4State_Init &&
	   fStateManager->GetCurrentState() != G4State_Idle));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

