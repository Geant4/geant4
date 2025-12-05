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
// gpaterno, October 2025
//
/// \file DetectorConstructionMessenger.hh
/// \brief Description of the DetectorConstruction messenger class
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
    
    G4UIcmdWithABool* fHybridSourceCmd{nullptr};
    
    G4UIdirectory* fCmdDir{nullptr};   
    G4UIcmdWithAString* fCrystalMaterialCmd{nullptr};
    G4UIcmdWith3VectorAndUnit* fCrystalSizeCmd{nullptr};
    G4UIcmdWithAString* fCrystalLatticeCmd{nullptr};
    G4UIcmdWithADouble* fCrystalAngleXCmd{nullptr};
    G4UIcmdWithADouble* fCrystalAngleYCmd{nullptr};
    G4UIcmdWithADouble* fCrystalBendingAngleCmd{nullptr};
    G4UIcmdWithABool* fRadModelCmd{nullptr};
    G4UIcmdWithABool* fOCeffectsCmd{nullptr};
    G4UIcmdWithAString* fPotentialPathCmd{nullptr};
    
    G4UIcmdWithADoubleAndUnit* fRadiatorConverterSepDistanceCmd{nullptr};
    G4UIcmdWith3VectorAndUnit* fConverterSizeCmd{nullptr};
    G4UIcmdWithABool* fGranularConverterCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fSphereRadiusCmd{nullptr};
    G4UIcmdWithAString* fConverterMaterialCmd{nullptr};
    
    G4UIcmdWithABool* fMagneticFieldCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fFieldValueCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fFieldRegionLengthCmd{nullptr};
    
    G4UIcmdWithABool* fCollimatorCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fCollimatorApertureCmd{nullptr};
    G4UIcmdWithAString* fCollimatorHoleCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fCollimatorThicknessCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fCollimatorSideCmd{nullptr};
    G4UIcmdWithADoubleAndUnit* fRadiatorCollimatorSepDistanceCmd{nullptr};
    
    G4UIcmdWith3VectorAndUnit* fVirtualDetectorSizeCmd{nullptr}; 
    
    G4UIcmdWithABool* fScoringCrystalExitCmd{nullptr}; 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

