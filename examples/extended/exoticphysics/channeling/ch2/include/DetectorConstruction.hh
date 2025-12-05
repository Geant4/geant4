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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ios.hh"
#include "globals.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include <vector>

#include "G4Region.hh"
#include "G4PVPlacement.hh"

#include "G4ChannelingFastSimModel.hh"

#include "DetectorConstructionMessenger.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    //methods to set the Crystal features
    void SetCrystalMaterial(const G4String& val) {fCrystalMaterialStr = val;}
    void SetCrystalSize(G4ThreeVector val) {fCrystalSize = val;}
    void SetCrystalBendingAngle(G4double val) {fBendingAngle = val;}
    void SetCrystalLattice(const G4String& val) {fLattice = val;}
    void SetCrystalAngleX(G4double val) {fAngleX = val;}
    void SetCrystalAngleY(G4double val) {fAngleY = val;}
    void SetRadiationModel(G4bool val) {fActivateRadiationModel = val;}
    void SetChannelingModel(G4bool val) {fActivateChannelingModel = val;}

    void SetCrystallineUndulatorAmplitude(G4double val)
    {fCrystallineUndulatorAmplitude = val;}
    void SetCrystallineUndulatorPeriod(G4double val)
    {fCrystallineUndulatorPeriod = val;}
    void SetCrystallineUndulatorPhase(G4double val)
    {fCrystallineUndulatorPhase = val;}

    void SetPotentialPath(const G4String& path){fPotentialPath = path;}
    void SetCrystalInternalGeometryPath(const G4String& path){fCrystalInternalGeometryPath = path;}
    void SetVirtualCollimatorHalfSize(G4double val) {fVirtualCollimatorHalfSize = val;}
    void SetMinPhotonEnergy(G4double val) {fMinPhotonEnergy = val;}
    void SetMaxBKPhotonEnergyInSpectrum(G4double val) {fMaxPhotonEnergySpectrum = val;}
    void SetNBinsSpectrum(G4int val) {fNBinsSpectrum = val;}
    void SetSamplingPhotonsNumber(G4int val) {fSamplingPhotonsNumber = val;}
    void SetNSmallTrajectorySteps(G4int val) {fNSmallTrajectorySteps = val;}
    void SetRadiationAngleFactor(G4double val) {fRadiationAngleFactor = val;}

    void SetMinPhotonEnergyAddStat(G4double val) {fMinPhotonEnergyAddStat = val;}
    void SetMaxPhotonEnergyAddStat(G4double val) {fMaxPhotonEnergyAddStat = val;}
    void SetMultiplePhotonStatistics(G4int val) {fTimesPhotonStatistics = val;}

    void SetDetectorSize(G4ThreeVector val) {fDetectorSize = val;}
    void SetDetectorFrontPositionZ(G4double val) {fDetectorFrontPosZ = val;}

    void SetParticleMinKinEnergy(G4double val)   {fParticleMinKinEnergy = val;}
    void SetProtonMinKinEnergy(G4double val)     {fProtonMinKinEnergy = val;}
    void SetAntiprotonMinKinEnergy(G4double val) {fAntiprotonMinKinEnergy = val;}
    void SetPiPlusMinKinEnergy(G4double val)     {fPiPlusMinKinEnergy = val;}
    void SetPiMinusMinKinEnergy(G4double val)    {fPiMinusMinKinEnergy = val;}
    void SetElectronMinKinEnergy(G4double val)   {fElectronMinKinEnergy = val;}
    void SetPositronMinKinEnergy(G4double val)   {fPositronMinKinEnergy = val;}
    void SetMuPlusMinKinEnergy(G4double val)     {fMuPlusMinKinEnergy = val;}
    void SetMuMinusMinKinEnergy(G4double val)    {fMuMinusMinKinEnergy = val;}

    void SetLindhardAngles(G4double val)           {fLindhardAngles = val;}
    void SetLindhardAnglesProton(G4double val)     {fLindhardAnglesProton = val;}
    void SetLindhardAnglesAntiproton(G4double val) {fLindhardAnglesAntiproton = val;}
    void SetLindhardAnglesPiPlus(G4double val)     {fLindhardAnglesPiPlus = val;}
    void SetLindhardAnglesPiMinus(G4double val)    {fLindhardAnglesPiMinus = val;}
    void SetLindhardAnglesElectron(G4double val)   {fLindhardAnglesElectron = val;}
    void SetLindhardAnglesPositron(G4double val)   {fLindhardAnglesPositron = val;}
    void SetLindhardAnglesMuPlus(G4double val)     {fLindhardAnglesMuPlus = val;}
    void SetLindhardAnglesMuMinus(G4double val)    {fLindhardAnglesMuMinus = val;}

  private:
    DetectorConstructionMessenger* fMessenger;

      //crystal features
      G4LogicalVolume* fLogicCrystal{nullptr};
      G4String fCrystalMaterialStr = "G4_Si";
      G4Material* fCrystalMaterial{nullptr};
      G4ThreeVector fCrystalSize;
      G4double fBendingAngle = 0.;
      G4String fLattice;
      G4double fAngleX = 0.;
      G4double fAngleY = 0.;
      G4bool fActivateRadiationModel = false;
      G4bool fActivateChannelingModel = true;

      //Crystal undulator parameters; default 0 => no undulator
      G4double fCrystallineUndulatorAmplitude = 0.;
      G4double fCrystallineUndulatorPeriod = 0.;
      G4double fCrystallineUndulatorPhase = 0.;

      G4ThreeVector fDetectorSize;
      G4double fDetectorFrontPosZ = 1*CLHEP::m;

      G4String fPotentialPath = "";
      G4String fCrystalInternalGeometryPath = "";

      G4double fMinPhotonEnergy = 1*CLHEP::MeV; //G4BaierKatkov default value
      G4double fMaxPhotonEnergySpectrum = 1*CLHEP::GeV; //G4BaierKatkov default value
      G4int fNBinsSpectrum = 110; //G4BaierKatkov default value
      G4int fSamplingPhotonsNumber = 150; //G4BaierKatkov default value
      G4int fNSmallTrajectorySteps = 10000; //G4BaierKatkov default value
      G4double fRadiationAngleFactor = 4; //G4BaierKatkov default value

      G4double fVirtualCollimatorHalfSize = 10.; // infinite collimator size
      G4double fMinPhotonEnergyAddStat = 1*CLHEP::MeV;
      G4double fMaxPhotonEnergyAddStat = 20*CLHEP::MeV;
      G4int fTimesPhotonStatistics = 1;

      G4double fParticleMinKinEnergy   = 200.*CLHEP::MeV;//G4ChannelingFastSimModel default value
      G4double fProtonMinKinEnergy     = 200.*CLHEP::MeV;
      G4double fAntiprotonMinKinEnergy = 200.*CLHEP::MeV;
      G4double fPiPlusMinKinEnergy     = 200.*CLHEP::MeV;
      G4double fPiMinusMinKinEnergy    = 200.*CLHEP::MeV;
      G4double fElectronMinKinEnergy   = 200.*CLHEP::MeV;
      G4double fPositronMinKinEnergy   = 200.*CLHEP::MeV;
      G4double fMuPlusMinKinEnergy     = 200.*CLHEP::MeV;
      G4double fMuMinusMinKinEnergy    = 200.*CLHEP::MeV;

      G4double fLindhardAngles = 100; //G4ChannelingFastSimModel default value
      G4double fLindhardAnglesProton     = 100;
      G4double fLindhardAnglesAntiproton = 100;
      G4double fLindhardAnglesPiPlus     = 100;
      G4double fLindhardAnglesPiMinus    = 100;
      G4double fLindhardAnglesElectron   = 100;
      G4double fLindhardAnglesPositron   = 100;
      G4double fLindhardAnglesMuPlus     = 100;
      G4double fLindhardAnglesMuMinus    = 100;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
