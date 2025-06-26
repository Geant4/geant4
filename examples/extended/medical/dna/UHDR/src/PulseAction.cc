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
// author: hoang tran

#include "PulseAction.hh"

#include "PulseActionMessenger.hh"

#include "G4RunManager.hh"
#include "G4Scheduler.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

#include <memory>

G4Mutex PulseAction::gUHDRMutex = G4MUTEX_INITIALIZER;  // Le Tuan Anh: for autolock in MT mode
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PulseAction::PulseAction(const G4String& pulse, G4bool useHisto)
  : G4UserTrackingAction(), fFileName(pulse), fUseHistoInput(useHisto)
{
  fpPulseInfo = std::make_unique<PulseInfo>(0);
  fpMessenger = std::make_unique<PulseActionMessenger>(this);
  if (fUseHistoInput) {
    InitializeForHistoInput();
  }
  else {
    Initialize();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PulseInfo::PulseInfo(G4double delayedTime) : G4VUserPulseInfo(), fDelayedTime(delayedTime) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PulseInfo::GetDelayedTime() const
{
  return fDelayedTime;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PulseInfo::~PulseInfo() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PulseAction::~PulseAction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PulseAction::PreUserTrackingAction(const G4Track* pTrack)
{
  if (!fActivePulse || fFileName == "") {
    return;
  }
  if (fActivePulse && fFileName == "") {
    return;
  }
  if (fFileName != "" && !fActivePulse) {
    return;
  }

  if (pTrack->GetParentID() == 0) {
    fDelayedTime = RandomizeInPulse();
    fpPulseInfo = std::make_unique<PulseInfo>(fDelayedTime);
    if (fVerbose > 1) {
      G4cout << "Particle comes at : " << G4BestUnit(fpPulseInfo->GetDelayedTime(), "Time")
             << G4endl;
    }
    if (fLonggestDelayedTime < fDelayedTime) {
      fLonggestDelayedTime = fDelayedTime;
    }
  }
  auto pPulseInfo = new PulseInfo(*fpPulseInfo);
  ((G4Track*)pTrack)->SetUserInformation(pPulseInfo);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PulseAction::Interpolate(const std::array<G4double, 5>& data)
{
  G4double e1 = data[0];
  G4double e2 = data[1];
  G4double e = data[2];
  G4double xs1 = data[3];
  G4double xs2 = data[4];
  G4double value = 0.;

  if ((e2 - e1) != 0) {
    G4double d1 = xs1;
    G4double d2 = xs2;
    value = (d1 + (d2 - d1) * (e - e1) / (e2 - e1));
  }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PulseAction::Initialize()
{
  if (fFileName.empty()) {
    return;
  }

  fPulseVector = {0.};  // Commence avec zéro comme première valeur
  fPulseData.clear();
  fPulseLarger = 0;
  G4MUTEXLOCK(&gUHDRMutex);  // Le Tuan Anh: for autolock in MT mode

  std::ifstream inputFile(fFileName);
  if (!inputFile) {
    G4ExceptionDescription exception;
    exception << "pulse Shape file not found. Please, provide : " << fFileName;
    G4Exception("PulseAction::Initialize()", "PulseAction01", FatalException, exception);
  }

  G4double time, pulseAmplitude;
  while (inputFile >> time) {
    time *= CLHEP::us;  // Conversion en unités

    if (!(inputFile >> pulseAmplitude)) {
      break;
    }
    if (time != fPulseVector.back()) {
      fPulseVector.push_back(time);
    }
    fPulseData[time] = pulseAmplitude;
    fPulseLarger = std::max(fPulseLarger, time);
  }
  G4MUTEXUNLOCK(&gUHDRMutex);  // Le Tuan Anh: for autolock in MT mode
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PulseAction::InitializeForHistoInput()
{
  // L.T. Anh: read histogram for Randomizing pulse using CLHEP::RandGeneral
  std::ostringstream FileName;
  FileName << fFileName;
  if (fFileName == "") {
    return;
  }
  fPulseVector.clear();
  G4int nbins = 0;
  G4double Tmin = 0, Tmax = 0, pdata;
  G4MUTEXLOCK(&gUHDRMutex);
  std::ifstream input(FileName.str().c_str());
  if (!input.is_open()) {
    G4ExceptionDescription exception;
    exception << "pulse Shape file " << fFileName
              << " not found. Please, provide correct file name!!!";
    G4Exception("PulseAction::InitializeForHistoInput()", "PulseAction01", FatalException,
                exception);
  }

  G4String aline;
  while (std::getline(input, aline)) {
    if (aline.empty()) continue;
    std::istringstream issLine(aline);
    G4String firstWord;
    issLine >> firstWord;
    if (firstWord == "#") continue;
    if (firstWord == "nbins") {
      issLine >> nbins;
    }
    else if (firstWord == "Tmin") {
      G4String units = "";
      issLine >> Tmin >> units;
      fTmin = Tmin * G4UIcommand::ValueOf(units);
    }
    else if (firstWord == "Tmax") {
      G4String units = "";
      issLine >> Tmax >> units;
      fTmax = Tmax * G4UIcommand::ValueOf(units);
    }
    else {
      pdata = std::stod(firstWord);
      fPulseVector.push_back(pdata);
    }
  }
  input.close();
  G4MUTEXUNLOCK(&gUHDRMutex);
  if (nbins != (G4int)fPulseVector.size() || nbins == 0) {
    G4ExceptionDescription exception;
    exception << "Nbins =  " << nbins;
    if (nbins != 0) exception << " not equal to data-size = " << fPulseVector.size();
    exception << "!!! \nPlease check the content/format of input file.";
    G4Exception("PulseAction::InitializeForHistoInput()", "PulseAction01", FatalException,
                exception);
  }
  fRandGeneral = std::make_unique<CLHEP::RandGeneral>(&fPulseVector.at(0), nbins);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PulseAction::RandomizeInPulse()
{
  if (fUseHistoInput) {
    // L.T. Anh: Randomizing pulse using CLHEP::RandGeneral
    G4double r = fRandGeneral->shoot();
    G4double tp = (fTmax - fTmin) * r + fTmin;
    return tp;
  }
  else {
    const G4double minTime = 0.;
    const G4double maxTime = fPulseLarger;  // us

    G4double MaximumPulse = 0.;
    G4int nSteps = 50;
    G4double value(minTime);

    for (G4int i = 0; i < nSteps; i++) {
      G4double PulseNumber = PulseSpectrum(value);
      if (PulseNumber >= MaximumPulse) {
        MaximumPulse = PulseNumber;
      }
      value += maxTime / nSteps;
    }

    G4double selectedPulse;
    do {
      selectedPulse = G4UniformRand() * (maxTime - minTime);
    } while (G4UniformRand() * MaximumPulse > PulseSpectrum(selectedPulse));

    return selectedPulse;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double PulseAction::PulseSpectrum(G4double time)
{
  G4double pulse = 0.;
  G4double valueT1, valueT2, xs1, xs2;
  auto t2 = std::upper_bound(fPulseVector.begin(), fPulseVector.end(), time);
  auto t1 = t2 - 1;
  valueT1 = *t1;
  valueT2 = *t2;
  xs1 = fPulseData[valueT1];
  xs2 = fPulseData[valueT2];
  auto xsProduct = xs1 * xs2;
  if (xsProduct != 0.) {
    std::array<G4double, 5> a = {valueT1, valueT2, time, xs1, xs2};
    pulse = Interpolate(a);
  }
  return pulse;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PulseAction::GetLonggestDelayedTime() const
{
  return fLonggestDelayedTime;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double PulseAction::GetPulseLarger() const
{
  return fPulseLarger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......