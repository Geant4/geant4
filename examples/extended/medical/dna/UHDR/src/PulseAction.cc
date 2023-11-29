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

#include <memory>

#include "PulseAction.hh"
#include "G4Track.hh"
#include "Randomize.hh"
#include "PulseActionMessenger.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PulseAction::PulseAction() :
    G4UserTrackingAction() {
  fpPulseInfo = std::make_unique<PulseInfo>(0);
  fpMessenger = std::make_unique<PulseActionMessenger>(this);
  Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PulseInfo::PulseInfo(G4double delayedTime)
    : G4VUserPulseInfo(), fDelayedTime(delayedTime) {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PulseInfo::GetDelayedTime() const {
  return fDelayedTime;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PulseInfo::~PulseInfo() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PulseAction::~PulseAction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PulseAction::PreUserTrackingAction(const G4Track *pTrack) {
  if(fActivePulse)
  {
    if (pTrack->GetParentID() == 0) {
      fDelayedTime = RandomizeInPulse();
      fpPulseInfo = std::make_unique<PulseInfo>(fDelayedTime);

      G4cout<<"Particle comes at : "<<G4BestUnit(fpPulseInfo->GetDelayedTime(),"Time")<<G4endl;
      if (fLonggestDelayedTime < fDelayedTime) {
        fLonggestDelayedTime = fDelayedTime;
      }
    }
    auto pPulseInfo = new PulseInfo(*fpPulseInfo);
    ((G4Track *) pTrack)->SetUserInformation(pPulseInfo);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PulseAction::Interpolate(const std::array<G4double,5>& data){
  G4double e1 = data[0];
  G4double e2 = data[1];
  G4double e = data[2];
  G4double xs1 = data[3];
  G4double xs2 = data[4];
  G4double value = 0.;
  if ((std::log10(e2) - std::log10(e1)) != 0) {
    G4double a = (std::log10(xs2) - std::log10(xs1))
                 / (std::log10(e2) - std::log10(e1));
    G4double b = std::log10(xs2) - a * std::log10(e2);
    G4double sigma = a * std::log10(e) + b;
    value = (std::pow(10., sigma));
  }

  if ((e2 - e1) != 0) {
    G4double d1 = xs1;
    G4double d2 = xs2;
    value = (d1 + (d2 - d1) * (e - e1) / (e2 - e1));
  }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PulseAction::Initialize() {
  std::ostringstream FileName;
  FileName << "pulseShape.dat";
  std::ifstream input(FileName.str().c_str());

  if (!input.is_open()) {
    G4ExceptionDescription exception;
    exception << "pulseShape.dat file not found. Please, provide";
    G4Exception("PulseAction::Initialize()", "PulseAction01",
                FatalException, exception);
  }

  fPulseVector.clear();
  fPulseVector.push_back(0.);
  while (!input.eof()) {
    double aTDummy;
    double pTDummy;
    input >> aTDummy;
    if (aTDummy != fPulseVector.back()) {
      fPulseVector.push_back(aTDummy);
    }
    input >> pTDummy;
    fPulseData[aTDummy] = pTDummy;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double PulseAction::RandomizeInPulse() {
  const G4double minTime = 0.;
  const G4double maxTime = fPulseLarger;// ns

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

  G4double selectedPulse = 0.;
  do {
    selectedPulse = G4UniformRand() * (maxTime - minTime);
  } while (G4UniformRand() * MaximumPulse > PulseSpectrum(selectedPulse));

  return selectedPulse;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double PulseAction::PulseSpectrum(G4double time) {
  G4double pulse = 0.;
  G4double valueT1 = 0;
  G4double valueT2 = 0;
  G4double xs1 = 0;
  G4double xs2 = 0;
  auto t2 = std::upper_bound(fPulseVector.begin(), fPulseVector.end(), time);
  auto t1 = t2 - 1;
  valueT1 = *t1;
  valueT2 = *t2;
  xs1 = fPulseData[valueT1];
  xs2 = fPulseData[valueT2];
  G4double xsProduct = xs1 * xs2;
  if (xsProduct != 0.) {
    std::array<G4double, 5> a = {valueT1, valueT2, time,xs1,xs2};
    pulse = Interpolate(a);
  }
  return pulse;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PulseAction::GetLonggestDelayedTime() const {
  return fLonggestDelayedTime;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
