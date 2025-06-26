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

#ifndef PULSE_PULSEACTION_HH
#define PULSE_PULSEACTION_HH 1

#include "CLHEP/Random/RandGeneral.h"

#include "G4AutoLock.hh"
#include "G4UserTrackingAction.hh"
#include "G4VUserPulseInfo.hh"
#include "G4VUserTrackInformation.hh"

#include <map>
#include <vector>

class G4ParticleDefinition;

class PulseActionMessenger;

class PulseInfo : public G4VUserPulseInfo
{
  public:
    explicit PulseInfo(G4double delayedTime);

    ~PulseInfo() override;

    G4double GetDelayedTime() const override;

  private:
    G4double fDelayedTime = 0;
};

class PulseAction : public G4UserTrackingAction
{
  public:
    using PulseMap = std::map<G4double, G4double>;

    PulseAction(const G4String& pulse, G4bool useHisto = false);

    ~PulseAction() override;

    void Initialize();

    void PreUserTrackingAction(const G4Track*) override;

    G4double RandomizeInPulse();

    G4double PulseSpectrum(G4double);

    static G4double Interpolate(const std::array<G4double, 5>& dat);

    G4double GetLonggestDelayedTime() const;

    inline void SetPulse(const G4bool& pulse) { fActivePulse = pulse; }

    inline G4bool IsActivedPulse() const { return fActivePulse; }
    // L.T. Anh added getter/setter for interpulse class:
    void SetLonggestDelayedTime(G4double lt) { fLonggestDelayedTime = lt; }
    const G4String GetPulseFileName() { return fFileName; }
    G4int GetVerbose() { return fVerbose; }
    G4double GetPulseLarger() const;

  protected:
    G4int fVerbose = 1;
    G4double fDelayedTime = 0;

private:
    void InitializeForHistoInput();
    std::unique_ptr<PulseInfo> fpPulseInfo;
    G4double fPulseLarger = 0;
    PulseMap fPulseData;
    std::vector<G4double> fPulseVector;
    G4double fLonggestDelayedTime = 0;
    std::unique_ptr<PulseActionMessenger> fpMessenger;
    G4bool fActivePulse = false;
    G4String fFileName = "";
    std::unique_ptr<CLHEP::RandGeneral> fRandGeneral{
      nullptr};  // L. T. Anh: pointer to radomize histogram from CLHEP
    G4bool fUseHistoInput{false};  // L. T. Anh: flag to use histo Input
    G4double fTmin{0.}, fTmax{0.};
    static G4Mutex gUHDRMutex;  // Le Tuan Anh: protect reading input files in MT mode
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif
