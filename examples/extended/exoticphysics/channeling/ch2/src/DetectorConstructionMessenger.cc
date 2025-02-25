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
/// \file DetectorConstructionMessenger.cc
/// \brief Implementation of the DetectorConstruction messenger class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstructionMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"

#include "G4RunManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstructionMessenger::DetectorConstructionMessenger(DetectorConstruction* det):
fDetector(det)
{
    fCmdDir = new G4UIdirectory("/crystal/");
    fCmdDir->SetGuidance("crystal Control");
    
    fCrystalMaterialCmd = new G4UIcmdWithAString("/crystal/setCrystalMaterial",this);
    fCrystalMaterialCmd->SetGuidance("Set Crystal Material");
    fCrystalMaterialCmd->SetParameterName("matname",false);
    fCrystalMaterialCmd->SetDefaultValue("G4_Si");
               
    fCrystalSizeCmd = new G4UIcmdWith3VectorAndUnit("/crystal/setCrystalSize",this);
    fCrystalSizeCmd->SetGuidance("Set Crystal size");
    fCrystalSizeCmd->SetParameterName("dimCrX","dimCrY","dimCrZ",false);
    fCrystalSizeCmd->SetUnitCategory("Length");
    fCrystalSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fCrystalLatticeCmd = new G4UIcmdWithAString("/crystal/setCrystalLattice",this);
    fCrystalLatticeCmd->
     SetGuidance("Set Crystal Lattice, use brackets (...) for planes and <...> for axes");
    fCrystalLatticeCmd->SetParameterName("lattice",false);
    fCrystalLatticeCmd->SetDefaultValue("(111)");
      
    fCrystalAngleXCmd = new G4UIcmdWithADoubleAndUnit("/crystal/setCrystalAngleX",this);
    fCrystalAngleXCmd->SetGuidance("Set crystal orientation with respect to the beam");
    fCrystalAngleXCmd->SetUnitCategory("Angle");
    fCrystalAngleXCmd->SetParameterName("angX",false);
    fCrystalAngleXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fCrystalAngleYCmd = new G4UIcmdWithADoubleAndUnit("/crystal/setCrystalAngleY",this);
    fCrystalAngleYCmd->SetGuidance("Set crystal orientation with respect to the beam");
    fCrystalAngleYCmd->SetUnitCategory("Angle");
    fCrystalAngleYCmd->SetParameterName("angY",false);
    fCrystalAngleYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fCrystalBendingAngleCmd =
            new G4UIcmdWithADoubleAndUnit("/crystal/setCrystalBendingAngle",this);
    fCrystalBendingAngleCmd->SetGuidance("Set crystal bending angle");
    fCrystalBendingAngleCmd->SetParameterName("bendingAngle",false);
    fCrystalBendingAngleCmd->SetUnitCategory("Angle");
    fCrystalBendingAngleCmd->SetRange("bendingAngle >= 0");
    fCrystalBendingAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fCrystallineUndulatorAmplitudeCmd =
            new G4UIcmdWithADoubleAndUnit
            ("/crystal/setCrystallineUndulatorAmplitude",this);
    fCrystallineUndulatorAmplitudeCmd->
            SetGuidance("Set crystalline undulator amplitude");
    fCrystallineUndulatorAmplitudeCmd->SetUnitCategory("Length");
    fCrystallineUndulatorAmplitudeCmd->
            SetParameterName("CrystallineUndulatorAmplitude",false);
    fCrystallineUndulatorAmplitudeCmd->
            SetRange("CrystallineUndulatorAmplitude > 0");
    fCrystallineUndulatorAmplitudeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fCrystallineUndulatorPeriodCmd =
            new G4UIcmdWithADoubleAndUnit("/crystal/setCrystallineUndulatorPeriod",this);
    fCrystallineUndulatorPeriodCmd->
            SetGuidance("Set crystalline undulator Period");
    fCrystallineUndulatorPeriodCmd->SetUnitCategory("Length");
    fCrystallineUndulatorPeriodCmd->
            SetParameterName("CrystallineUndulatorPeriod",false);
    fCrystallineUndulatorPeriodCmd->
            SetRange("CrystallineUndulatorPeriod > 0");
    fCrystallineUndulatorPeriodCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fCrystallineUndulatorPhaseCmd =
            new G4UIcmdWithADouble("/crystal/setCrystallineUndulatorPhase",this);
    fCrystallineUndulatorPhaseCmd->
            SetGuidance("Set crystalline undulator phase");
    fCrystallineUndulatorPhaseCmd->
            SetParameterName("CrystallineUndulatorPhase",false);
    fCrystallineUndulatorPhaseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fDetectorSizeCmd = new G4UIcmdWith3VectorAndUnit("/crystal/setDetectorSize",this);
    fDetectorSizeCmd->SetGuidance("Set detector size");
    fDetectorSizeCmd->SetParameterName("dimDetX","dimDetY","dimDetZ",false);
    fDetectorSizeCmd->SetUnitCategory("Length");
    fDetectorSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fDetectorFrontPosZCmd =
            new G4UIcmdWithADoubleAndUnit("/crystal/setFrontPositionZ",this);
    fDetectorFrontPosZCmd->SetGuidance("Set detector front position Z");
    fDetectorFrontPosZCmd->SetParameterName("frontPosDetZ",false);
    fDetectorFrontPosZCmd->SetUnitCategory("Length");
    fDetectorFrontPosZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fPotentialPathCmd = new G4UIcmdWithAString("/crystal/setChannelingDataPath",this);
    fPotentialPathCmd->
            SetGuidance("Set the path where to find the available data "
                        "for the G4ChannelingFastSimModel "
                        "if different from G4CHANNELINGDATA");
    fPotentialPathCmd->SetParameterName("channelingDataPath",false);
    fPotentialPathCmd->SetDefaultValue("");

    fChannelingModelCmd = new G4UIcmdWithABool("/crystal/setChannelingModel", this);
    fChannelingModelCmd->SetGuidance("Activate/deactivate G4ChannelingFastSimModel");
    fChannelingModelCmd->SetParameterName("ChannelingModel",true);
    fChannelingModelCmd->SetDefaultValue(false);

    fRadModelCmd = new G4UIcmdWithABool("/crystal/setRadiationModel", this);
    fRadModelCmd->SetGuidance("Activate/deactivate G4BaierKatkov");
    fRadModelCmd->SetParameterName("ActivateRadiationModel",true);
    fRadModelCmd->SetDefaultValue(false);

    fMinPhotonEnergyCmd =
            new G4UIcmdWithADoubleAndUnit("/crystal/setMinPhotonEnergy",this);
    fMinPhotonEnergyCmd->
            SetGuidance("Set the low energy threshold for "
                        "the photons emitted in G4BaierKatkov");
    fMinPhotonEnergyCmd->SetParameterName("MinPhotonEnergy",false);
    fMinPhotonEnergyCmd->SetUnitCategory("Energy");
    fMinPhotonEnergyCmd->SetRange("MinPhotonEnergy > 0");
    fMinPhotonEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fSamplingPhotonsNumberCmd =
            new G4UIcmdWithAnInteger("/crystal/setSamplingPhotonsNumber",this);
    fSamplingPhotonsNumberCmd->
            SetGuidance("Set SamplingPhotonsNumber in G4BaierKatkov");
    fSamplingPhotonsNumberCmd->SetParameterName("SamplingPhotonsNumber",false);
    fSamplingPhotonsNumberCmd->SetRange("SamplingPhotonsNumber>1");
    fSamplingPhotonsNumberCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fNSmallTrajectoryStepsCmd =
            new G4UIcmdWithAnInteger("/crystal/setNSmallTrajectorySteps",this);
    fNSmallTrajectoryStepsCmd->
            SetGuidance("Set NSmallTrajectorySteps in G4BaierKatkov");
    fNSmallTrajectoryStepsCmd->SetParameterName("NSmallTrajectorySteps",false);
    fNSmallTrajectoryStepsCmd->SetRange("NSmallTrajectorySteps>1");
    fNSmallTrajectoryStepsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fRadiationAngleFactorCmd =
            new G4UIcmdWithADouble("/crystal/setRadiationAngleFactor",this);
    fRadiationAngleFactorCmd->SetGuidance("Set Radiation Angle Factor");
    fRadiationAngleFactorCmd->SetParameterName("RadiationAngleFactor",false);
    fRadiationAngleFactorCmd->SetRange("RadiationAngleFactor > 0");
    fRadiationAngleFactorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fMinPhotonEnergyAddStatCmd =
    new G4UIcmdWithADoubleAndUnit("/crystal/AddPhotonStatistics/setMinPhotonEnergy",this);
    fMinPhotonEnergyAddStatCmd->
            SetGuidance("Set the min energy in the range to increase "
                        "the sampling photon statistics in G4BaierKatkov");
    fMinPhotonEnergyAddStatCmd->SetParameterName("addStatMinPhotonEnergy",false);
    fMinPhotonEnergyAddStatCmd->SetUnitCategory("Energy");
    fMinPhotonEnergyAddStatCmd->SetRange("addStatMinPhotonEnergy > 0");
    fMinPhotonEnergyAddStatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fMaxPhotonEnergyAddStatCmd =
    new G4UIcmdWithADoubleAndUnit("/crystal/AddPhotonStatistics/setMaxPhotonEnergy",this);
    fMaxPhotonEnergyAddStatCmd->
            SetGuidance("Set the max energy in the range to increase "
                        "the sampling photon statistics in G4BaierKatkov");
    fMaxPhotonEnergyAddStatCmd->SetParameterName("addStatMaxPhotonEnergy",false);
    fMaxPhotonEnergyAddStatCmd->SetUnitCategory("Energy");
    fMaxPhotonEnergyAddStatCmd->SetRange("addStatMaxPhotonEnergy > 0");
    fMaxPhotonEnergyAddStatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fTimesPhotonStatisticsCmd =
    new G4UIcmdWithAnInteger("/crystal/AddPhotonStatistics/setMultiplePhotonStatistics",
                             this);
    fTimesPhotonStatisticsCmd->
           SetGuidance("Set multiple of the sampling photon statistics in G4BaierKatkov");
    fTimesPhotonStatisticsCmd->SetParameterName("timesPhotonStatistics",false);
    fTimesPhotonStatisticsCmd->SetRange("timesPhotonStatistics > 1");
    fTimesPhotonStatisticsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fParticleMinKinEnergyCmd =
            new G4UIcmdWithADoubleAndUnit("/crystal/setParticleMinKinEnergy",this);
    fParticleMinKinEnergyCmd->
            SetGuidance("Set the low energy threshold for particle "
                        "to enter the G4ChannelingFastSimModel");
    fParticleMinKinEnergyCmd->SetParameterName("partLEth",false);
    fParticleMinKinEnergyCmd->SetUnitCategory("Energy");
    fParticleMinKinEnergyCmd->SetRange("partLEth > 0");
    fParticleMinKinEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fProtonMinKinEnergyCmd =
            new G4UIcmdWithADoubleAndUnit("/crystal/setParticleMinKinEnergy/proton",this);
    fProtonMinKinEnergyCmd->
            SetGuidance("Set the low energy threshold for proton "
                        "to enter the G4ChannelingFastSimModel");
    fProtonMinKinEnergyCmd->SetParameterName("protonLEth",false);
    fProtonMinKinEnergyCmd->SetUnitCategory("Energy");
    fProtonMinKinEnergyCmd->SetRange("protonLEth > 0");
    fProtonMinKinEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAntiprotonMinKinEnergyCmd =
       new G4UIcmdWithADoubleAndUnit("/crystal/setParticleMinKinEnergy/anti_proton",this);
    fAntiprotonMinKinEnergyCmd->SetGuidance("Set the low energy threshold for anti_proton"
                                            " to enter the G4ChannelingFastSimModel");
    fAntiprotonMinKinEnergyCmd->SetParameterName("anti_protonLEth",false);
    fAntiprotonMinKinEnergyCmd->SetUnitCategory("Energy");
    fAntiprotonMinKinEnergyCmd->SetRange("anti_protonLEth > 0");
    fAntiprotonMinKinEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fPiPlusMinKinEnergyCmd =
            new G4UIcmdWithADoubleAndUnit("/crystal/setParticleMinKinEnergy/pi+",this);
    fPiPlusMinKinEnergyCmd->
            SetGuidance("Set the low energy threshold for pi+ "
                        "to enter the G4ChannelingFastSimModel");
    fPiPlusMinKinEnergyCmd->SetParameterName("piPlusLEth",false);
    fPiPlusMinKinEnergyCmd->SetUnitCategory("Energy");
    fPiPlusMinKinEnergyCmd->SetRange("piPlusLEth > 0");
    fPiPlusMinKinEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fPiMinusMinKinEnergyCmd =
            new G4UIcmdWithADoubleAndUnit("/crystal/setParticleMinKinEnergy/pi-",this);
    fPiMinusMinKinEnergyCmd->SetGuidance("Set the low energy threshold for pi- "
                                         "to enter the G4ChannelingFastSimModel");
    fPiMinusMinKinEnergyCmd->SetParameterName("piMinusLEth",false);
    fPiMinusMinKinEnergyCmd->SetUnitCategory("Energy");
    fPiMinusMinKinEnergyCmd->SetRange("piMinusLEth > 0");
    fPiMinusMinKinEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fPositronMinKinEnergyCmd =
            new G4UIcmdWithADoubleAndUnit("/crystal/setParticleMinKinEnergy/e+",this);
    fPositronMinKinEnergyCmd->SetGuidance("Set the low energy threshold for e+ "
                                          "to enter the G4ChannelingFastSimModel");
    fPositronMinKinEnergyCmd->SetParameterName("ePlusLEth",false);
    fPositronMinKinEnergyCmd->SetUnitCategory("Energy");
    fPositronMinKinEnergyCmd->SetRange("ePlusLEth > 0");
    fPositronMinKinEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fElectronMinKinEnergyCmd =
            new G4UIcmdWithADoubleAndUnit("/crystal/setParticleMinKinEnergy/e-",this);
    fElectronMinKinEnergyCmd->SetGuidance("Set the low energy threshold for e- "
                                          "to enter the G4ChannelingFastSimModel");
    fElectronMinKinEnergyCmd->SetParameterName("eMinusLEth",false);
    fElectronMinKinEnergyCmd->SetUnitCategory("Energy");
    fElectronMinKinEnergyCmd->SetRange("eMinusLEth > 0");
    fElectronMinKinEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fMuPlusMinKinEnergyCmd =
            new G4UIcmdWithADoubleAndUnit("/crystal/setParticleMinKinEnergy/mu+",this);
    fMuPlusMinKinEnergyCmd->SetGuidance("Set the low energy threshold for mu+ "
                                        "to enter the G4ChannelingFastSimModel");
    fMuPlusMinKinEnergyCmd->SetParameterName("muPlusLEth",false);
    fMuPlusMinKinEnergyCmd->SetUnitCategory("Energy");
    fMuPlusMinKinEnergyCmd->SetRange("muPlusLEth > 0");
    fMuPlusMinKinEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fMuMinusMinKinEnergyCmd =
            new G4UIcmdWithADoubleAndUnit("/crystal/setParticleMinKinEnergy/mu-",this);
    fMuMinusMinKinEnergyCmd->SetGuidance("Set the low energy threshold for mu- "
                                         "to enter the G4ChannelingFastSimModel");
    fMuMinusMinKinEnergyCmd->SetParameterName("muMinusLEth",false);
    fMuMinusMinKinEnergyCmd->SetUnitCategory("Energy");
    fMuMinusMinKinEnergyCmd->SetRange("muMinusLEth > 0");
    fMuMinusMinKinEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    fLindhardAnglesCmd = new G4UIcmdWithADouble("/crystal/setLindhardAngles",this);
    fLindhardAnglesCmd->
            SetGuidance("Set high angular threshold for particle to enter "
                        "the G4ChannelingFastSimModel expressed in Lindhard angles");
    fLindhardAnglesCmd->SetParameterName("LindhardAngles",false);
    fLindhardAnglesCmd->SetRange("LindhardAngles >= 0");
    fLindhardAnglesCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fLindhardAnglesProtonCmd =
            new G4UIcmdWithADouble("/crystal/setLindhardAngles/proton",this);
    fLindhardAnglesProtonCmd->
            SetGuidance("Set high angular threshold for proton to enter "
                        "the G4ChannelingFastSimModel expressed in Lindhard angles");
    fLindhardAnglesProtonCmd->SetParameterName("LindhardAnglesProton",false);
    fLindhardAnglesProtonCmd->SetRange("LindhardAnglesProton >= 0");
    fLindhardAnglesProtonCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fLindhardAnglesAntiprotonCmd =
            new G4UIcmdWithADouble("/crystal/setLindhardAngles/anti_proton",this);
    fLindhardAnglesAntiprotonCmd->
            SetGuidance("Set high angular threshold for anti_proton to enter "
                        "the G4ChannelingFastSimModel expressed in Lindhard angles");
    fLindhardAnglesAntiprotonCmd->SetParameterName("LindhardAnglesAnti_proton",false);
    fLindhardAnglesAntiprotonCmd->SetRange("LindhardAnglesAnti_proton >= 0");
    fLindhardAnglesAntiprotonCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fLindhardAnglesPiPlusCmd =
            new G4UIcmdWithADouble("/crystal/setLindhardAngles/pi+",this);
    fLindhardAnglesPiPlusCmd->
            SetGuidance("Set high angular threshold for pi+ to enter "
                        "the G4ChannelingFastSimModel expressed in Lindhard angles");
    fLindhardAnglesPiPlusCmd->SetParameterName("LindhardAnglesPiPlus",false);
    fLindhardAnglesPiPlusCmd->SetRange("LindhardAnglesPiPlus >= 0");
    fLindhardAnglesPiPlusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fLindhardAnglesPiMinusCmd =
            new G4UIcmdWithADouble("/crystal/setLindhardAngles/pi-",this);
    fLindhardAnglesPiMinusCmd->
            SetGuidance("Set high angular threshold for pi- to enter "
                        "the G4ChannelingFastSimModel expressed in Lindhard angles");
    fLindhardAnglesPiMinusCmd->SetParameterName("LindhardAnglesPiMinus",false);
    fLindhardAnglesPiMinusCmd->SetRange("LindhardAnglesPiMinus >= 0");
    fLindhardAnglesPiMinusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fLindhardAnglesPositronCmd =
            new G4UIcmdWithADouble("/crystal/setLindhardAngles/e+",this);
    fLindhardAnglesPositronCmd->
            SetGuidance("Set high angular threshold for e+ to enter "
                        "the G4ChannelingFastSimModel expressed in Lindhard angles");
    fLindhardAnglesPositronCmd->SetParameterName("LindhardAnglesPositron",false);
    fLindhardAnglesPositronCmd->SetRange("LindhardAnglesPositron >= 0");
    fLindhardAnglesPositronCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fLindhardAnglesElectronCmd =
            new G4UIcmdWithADouble("/crystal/setLindhardAngles/e-",this);
    fLindhardAnglesElectronCmd->
            SetGuidance("Set high angular threshold for e- to enter "
                        "the G4ChannelingFastSimModel expressed in Lindhard angles");
    fLindhardAnglesElectronCmd->SetParameterName("LindhardAnglesElectron",false);
    fLindhardAnglesElectronCmd->SetRange("LindhardAnglesElectron >= 0");
    fLindhardAnglesElectronCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fLindhardAnglesMuPlusCmd =
            new G4UIcmdWithADouble("/crystal/setLindhardAngles/mu+",this);
    fLindhardAnglesMuPlusCmd->
            SetGuidance("Set high angular threshold for mu+ to enter "
                        "the G4ChannelingFastSimModel expressed in Lindhard angles");
    fLindhardAnglesMuPlusCmd->SetParameterName("LindhardAnglesMuPlus",false);
    fLindhardAnglesMuPlusCmd->SetRange("LindhardAnglesMuPlus >= 0");
    fLindhardAnglesMuPlusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fLindhardAnglesMuMinusCmd =
            new G4UIcmdWithADouble("/crystal/setLindhardAngles/mu-",this);
    fLindhardAnglesMuMinusCmd->
            SetGuidance("Set high angular threshold for mu- to enter "
                        "the G4ChannelingFastSimModel expressed in Lindhard angles");
    fLindhardAnglesMuMinusCmd->SetParameterName("LindhardAnglesMuMinus",false);
    fLindhardAnglesMuMinusCmd->SetRange("LindhardAnglesMuMinus >= 0");
    fLindhardAnglesMuMinusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstructionMessenger::~DetectorConstructionMessenger()
{
    delete fCmdDir;
    delete fCrystalMaterialCmd;
    delete fCrystalSizeCmd;
    delete fCrystalLatticeCmd;
    delete fCrystalAngleXCmd;
    delete fCrystalAngleYCmd;
    delete fCrystalBendingAngleCmd;
    delete fRadModelCmd;
    delete fChannelingModelCmd;
    
    delete fCrystallineUndulatorAmplitudeCmd;
    delete fCrystallineUndulatorPeriodCmd;
    delete fCrystallineUndulatorPhaseCmd;
    
    delete fDetectorSizeCmd;
    delete fDetectorFrontPosZCmd;

    delete fPotentialPathCmd;
    
    delete fMinPhotonEnergyCmd;
    delete fSamplingPhotonsNumberCmd;
    delete fNSmallTrajectoryStepsCmd;
    delete fRadiationAngleFactorCmd;
    delete fMinPhotonEnergyAddStatCmd;
    delete fMaxPhotonEnergyAddStatCmd;
    delete fTimesPhotonStatisticsCmd;

    delete fParticleMinKinEnergyCmd;
    delete fProtonMinKinEnergyCmd;
    delete fAntiprotonMinKinEnergyCmd;
    delete fPiPlusMinKinEnergyCmd;
    delete fPiMinusMinKinEnergyCmd;
    delete fElectronMinKinEnergyCmd;
    delete fPositronMinKinEnergyCmd;
    delete fMuPlusMinKinEnergyCmd;
    delete fMuMinusMinKinEnergyCmd;

    delete fLindhardAnglesCmd;
    delete fLindhardAnglesProtonCmd;
    delete fLindhardAnglesAntiprotonCmd;
    delete fLindhardAnglesPiPlusCmd;
    delete fLindhardAnglesPiMinusCmd;
    delete fLindhardAnglesElectronCmd;
    delete fLindhardAnglesPositronCmd;
    delete fLindhardAnglesMuPlusCmd;
    delete fLindhardAnglesMuMinusCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstructionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == fCrystalMaterialCmd)
        {fDetector->SetCrystalMaterial(newValue);}
    if (command == fCrystalSizeCmd) 
        {fDetector->SetCrystalSize(fCrystalSizeCmd->GetNew3VectorValue(newValue));}
    if (command == fCrystalLatticeCmd)
        {fDetector->SetCrystalLattice(newValue);}
    if (command == fCrystalAngleXCmd) 
        {fDetector->SetCrystalAngleX(fCrystalAngleXCmd->GetNewDoubleValue(newValue));}
    if (command == fCrystalAngleYCmd) 
        {fDetector->SetCrystalAngleY(fCrystalAngleYCmd->GetNewDoubleValue(newValue));}
    if (command == fCrystalBendingAngleCmd) 
        {fDetector->SetCrystalBendingAngle(
                    fCrystalBendingAngleCmd->GetNewDoubleValue(newValue));}
    if (command == fRadModelCmd) 
        {fDetector->SetRadiationModel(fRadModelCmd->GetNewBoolValue(newValue));}
    if (command == fChannelingModelCmd)
        {fDetector->SetChannelingModel(fChannelingModelCmd->GetNewBoolValue(newValue));}
    
    if (command == fCrystallineUndulatorAmplitudeCmd)
        {fDetector->SetCrystallineUndulatorAmplitude
                (fCrystallineUndulatorAmplitudeCmd->GetNewDoubleValue(newValue));}
    if (command == fCrystallineUndulatorPeriodCmd)
        {fDetector->SetCrystallineUndulatorPeriod
                (fCrystallineUndulatorPeriodCmd->GetNewDoubleValue(newValue));}
    if (command == fCrystallineUndulatorPhaseCmd)
        {fDetector->SetCrystallineUndulatorPhase
                (fCrystallineUndulatorPhaseCmd->GetNewDoubleValue(newValue));}

    if (command == fDetectorSizeCmd)
        {fDetector->SetDetectorSize(fDetectorSizeCmd->GetNew3VectorValue(newValue));}
    if (command == fDetectorFrontPosZCmd)
        {fDetector->SetDetectorFrontPositionZ(
                    fDetectorFrontPosZCmd->GetNewDoubleValue(newValue));}

    if (command == fPotentialPathCmd)
        {fDetector->SetPotentialPath(newValue);}
                
    if (command == fMinPhotonEnergyCmd)
        {fDetector->SetMinPhotonEnergy(
                    fMinPhotonEnergyCmd->GetNewDoubleValue(newValue));}
    if (command == fSamplingPhotonsNumberCmd)
        {fDetector->SetSamplingPhotonsNumber(
                    fSamplingPhotonsNumberCmd->GetNewIntValue(newValue));}
    if (command == fNSmallTrajectoryStepsCmd)
        {fDetector->SetNSmallTrajectorySteps(
                    fNSmallTrajectoryStepsCmd->GetNewIntValue(newValue));}
    if (command == fRadiationAngleFactorCmd)
        {fDetector->SetRadiationAngleFactor(
                    fRadiationAngleFactorCmd->GetNewDoubleValue(newValue));}

    if (command == fMinPhotonEnergyAddStatCmd)
        {fDetector->SetMinPhotonEnergyAddStat(
                    fMinPhotonEnergyAddStatCmd->GetNewDoubleValue(newValue));}
    if (command == fMaxPhotonEnergyAddStatCmd)
        {fDetector->SetMaxPhotonEnergyAddStat(
                    fMaxPhotonEnergyAddStatCmd->GetNewDoubleValue(newValue));}
    if (command == fTimesPhotonStatisticsCmd)
        {fDetector->SetMultiplePhotonStatistics(
                    fTimesPhotonStatisticsCmd->GetNewIntValue(newValue));}

    if (command == fParticleMinKinEnergyCmd)
        {fDetector->SetParticleMinKinEnergy(
                    fParticleMinKinEnergyCmd->GetNewDoubleValue(newValue));}
    if (command == fProtonMinKinEnergyCmd)
        {fDetector->SetProtonMinKinEnergy(
                    fProtonMinKinEnergyCmd->GetNewDoubleValue(newValue));}
    if (command == fAntiprotonMinKinEnergyCmd)
        {fDetector->SetAntiprotonMinKinEnergy(
                    fAntiprotonMinKinEnergyCmd->GetNewDoubleValue(newValue));}
    if (command == fPiPlusMinKinEnergyCmd)
        {fDetector->SetPiPlusMinKinEnergy(
                    fPiPlusMinKinEnergyCmd->GetNewDoubleValue(newValue));}
    if (command == fPiMinusMinKinEnergyCmd)
        {fDetector->SetPiMinusMinKinEnergy(
                    fPiMinusMinKinEnergyCmd->GetNewDoubleValue(newValue));}
    if (command == fElectronMinKinEnergyCmd)
        {fDetector->SetElectronMinKinEnergy(
                    fElectronMinKinEnergyCmd->GetNewDoubleValue(newValue));}
    if (command == fPositronMinKinEnergyCmd)
        {fDetector->SetPositronMinKinEnergy(
                    fPositronMinKinEnergyCmd->GetNewDoubleValue(newValue));}
    if (command == fMuPlusMinKinEnergyCmd)
        {fDetector->SetMuPlusMinKinEnergy(
                    fMuPlusMinKinEnergyCmd->GetNewDoubleValue(newValue));}
    if (command == fMuMinusMinKinEnergyCmd)
        {fDetector->SetMuMinusMinKinEnergy(
                    fMuMinusMinKinEnergyCmd->GetNewDoubleValue(newValue));}

    if (command == fLindhardAnglesCmd) 
        {fDetector->SetLindhardAngles(
                    fLindhardAnglesCmd->GetNewDoubleValue(newValue));}
    if (command == fLindhardAnglesProtonCmd)
        {fDetector->SetLindhardAnglesProton(
                    fLindhardAnglesProtonCmd->GetNewDoubleValue(newValue));}
    if (command == fLindhardAnglesAntiprotonCmd)
        {fDetector->SetLindhardAnglesAntiproton(
                    fLindhardAnglesAntiprotonCmd->GetNewDoubleValue(newValue));}
    if (command == fLindhardAnglesPiPlusCmd)
        {fDetector->SetLindhardAnglesPiPlus(
                    fLindhardAnglesPiPlusCmd->GetNewDoubleValue(newValue));}
    if (command == fLindhardAnglesPiMinusCmd)
        {fDetector->SetLindhardAnglesPiMinus(
                    fLindhardAnglesPiMinusCmd->GetNewDoubleValue(newValue));}
    if (command == fLindhardAnglesElectronCmd)
        {fDetector->SetLindhardAnglesElectron(
                    fLindhardAnglesElectronCmd->GetNewDoubleValue(newValue));}
    if (command == fLindhardAnglesPositronCmd)
        {fDetector->SetLindhardAnglesPositron(
                    fLindhardAnglesPositronCmd->GetNewDoubleValue(newValue));}
    if (command == fLindhardAnglesMuPlusCmd)
        {fDetector->SetLindhardAnglesMuPlus(
                    fLindhardAnglesMuPlusCmd->GetNewDoubleValue(newValue));}
    if (command == fLindhardAnglesMuMinusCmd)
        {fDetector->SetLindhardAnglesMuMinus(
                    fLindhardAnglesMuMinusCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
