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
/// \file DetectorConstructionMessenger.hh
/// \brief Description of the DetectorConstruction messenger class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstructionMessenger_h
#define DetectorConstructionMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Detector construction messenger class to define custom commands
/// to control the geometry and other settings.

class DetectorConstructionMessenger: public G4UImessenger
{
public:
    DetectorConstructionMessenger(DetectorConstruction* mpga);
    ~DetectorConstructionMessenger();

    void SetNewValue(G4UIcommand* command, G4String newValues) override;

private:
    DetectorConstruction* fDetector{nullptr};
    
    G4UIdirectory* fCmdDir{nullptr};   
    G4UIcmdWithAString* fCrystalMaterialCmd{nullptr};
    G4UIcmdWith3VectorAndUnit* fCrystalSizeCmd{nullptr};
    G4UIcmdWithAString* fCrystalLatticeCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fCrystalAngleXCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fCrystalAngleYCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fCrystalBendingAngleCmd{nullptr};
    G4UIcmdWithABool* fRadModelCmd{nullptr};
    G4UIcmdWithABool* fChannelingModelCmd{nullptr};
    
    G4UIcmdWith3VectorAndUnit* fDetectorSizeCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fDetectorFrontPosZCmd{nullptr};

    G4UIcmdWithADoubleAndUnit* fCrystallineUndulatorAmplitudeCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fCrystallineUndulatorPeriodCmd{nullptr};
    G4UIcmdWithADouble* fCrystallineUndulatorPhaseCmd{nullptr};
    
    G4UIcmdWithAString* fPotentialPathCmd{nullptr};
    
    G4UIcmdWithADoubleAndUnit* fMinPhotonEnergyCmd{nullptr};
    G4UIcmdWithAnInteger* fSamplingPhotonsNumberCmd{nullptr};
    G4UIcmdWithAnInteger* fNSmallTrajectoryStepsCmd{nullptr};
    G4UIcmdWithADouble* fRadiationAngleFactorCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fMinPhotonEnergyAddStatCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fMaxPhotonEnergyAddStatCmd{nullptr};
    G4UIcmdWithAnInteger* fTimesPhotonStatisticsCmd{nullptr};

    G4UIcmdWithADoubleAndUnit* fParticleMinKinEnergyCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fProtonMinKinEnergyCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fAntiprotonMinKinEnergyCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fPiPlusMinKinEnergyCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fPiMinusMinKinEnergyCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fElectronMinKinEnergyCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fPositronMinKinEnergyCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fMuPlusMinKinEnergyCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fMuMinusMinKinEnergyCmd{nullptr};

    G4UIcmdWithADouble* fLindhardAnglesCmd{nullptr};
    G4UIcmdWithADouble* fLindhardAnglesProtonCmd{nullptr};
    G4UIcmdWithADouble* fLindhardAnglesAntiprotonCmd{nullptr};
    G4UIcmdWithADouble* fLindhardAnglesPiPlusCmd{nullptr};
    G4UIcmdWithADouble* fLindhardAnglesPiMinusCmd{nullptr};
    G4UIcmdWithADouble* fLindhardAnglesElectronCmd{nullptr};
    G4UIcmdWithADouble* fLindhardAnglesPositronCmd{nullptr};
    G4UIcmdWithADouble* fLindhardAnglesMuPlusCmd{nullptr};
    G4UIcmdWithADouble* fLindhardAnglesMuMinusCmd{nullptr};

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

