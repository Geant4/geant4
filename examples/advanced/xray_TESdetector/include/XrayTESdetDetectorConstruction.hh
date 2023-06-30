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
/// \file XrayTESdetDetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
// Authors: P.Dondero (paolo.dondero@cern.ch), R.Stanzani (ronny.stanzani@cern.ch)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#ifndef XrayTESdetDetectorConstruction_h
#define XrayTESdetDetectorConstruction_h 1

class XrayTESdetDetectorMessenger;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4UserLimits;
class G4PVPlacement;

#include "G4VUserDetectorConstruction.hh"
#include "G4GDMLParser.hh"
#include "G4PVParameterised.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class XrayTESdetDetectorConstruction : public G4VUserDetectorConstruction
{
   public:

     explicit XrayTESdetDetectorConstruction();
     ~XrayTESdetDetectorConstruction() override;
     G4VPhysicalVolume* Construct() override;
     void SetMaxStep (G4double);
     void SetReadFile(G4String);

   private:

     // GDMLparser
     G4GDMLParser fParser;
     G4String fReadFile;
     G4UserLimits* fStepLimitPolyfilta;
     G4UserLimits* fStepLimitAlfilta;
     G4UserLimits* fStepLimitmembrane;
     G4PVPlacement* fExperimentalHall_phys;
     XrayTESdetDetectorMessenger* fDetectorMessenger;
     G4VPhysicalVolume* ConstructDetector();
     G4PVParameterised* fPhysiDet;

     // Logical volumes
     G4LogicalVolume* fExperimentalHall_log;
     G4LogicalVolume* fElement_log;
     G4LogicalVolume* fMembranepxl_log;
     G4LogicalVolume* fDetector_log;
     G4LogicalVolume* fBipxl_log;
     G4LogicalVolume* fACDpxl_log;
     G4LogicalVolume* fInACDpxl_log;
     G4LogicalVolume* fGridpiece_log;
     G4LogicalVolume* fWafer_log;
     G4LogicalVolume* fACplate_log;
     G4LogicalVolume* fCagewall_log;
     G4LogicalVolume* fMinicagewall_log;
     G4LogicalVolume* fExtcagewall_log;
     G4LogicalVolume* fExtboard_log;
     G4LogicalVolume* fMesh_log;
     G4LogicalVolume* fSphere_log;
     G4LogicalVolume* fBSC_log;
     G4LogicalVolume* fTrapezoid_A_log;
     G4LogicalVolume* fTrapezoid_B_log;
     G4LogicalVolume* fTrapezoid_C_log;
     G4LogicalVolume* fTrapezoid_D_log;
     G4LogicalVolume* fWorld_log;

     // Physical volumes
     G4VPhysicalVolume* fBipxl_phys;
     G4VPhysicalVolume* fMembranepxl_phys;
     G4VPhysicalVolume* fGridpiece_phys;
     G4VPhysicalVolume* fDetector_phys;
     G4VPhysicalVolume* fACDpxl_phys;
     G4VPhysicalVolume* fInACDpxl_phys;
     G4VPhysicalVolume* fWafer_phys;
     G4VPhysicalVolume* fACplate_phys;
     G4VPhysicalVolume* fCagecolumn_phys;
     G4VPhysicalVolume* fCagecolumn2_phys;
     G4VPhysicalVolume* fCagecolumn3_phys;
     G4VPhysicalVolume* fCagecolumn4_phys;
     G4VPhysicalVolume* fCagecolumn5_phys;
     G4VPhysicalVolume* fCagecolumn6_phys;
     G4VPhysicalVolume* fEmptycagecolumn_phys;
     G4VPhysicalVolume* fEmptycagecolumn2_phys;
     G4VPhysicalVolume* fEmptycagecolumn3_phys;
     G4VPhysicalVolume* fCagewall_phys;
     G4VPhysicalVolume* fCagewall2_phys;
     G4VPhysicalVolume* fCagewall3_phys;
     G4VPhysicalVolume* fCagewall4_phys;
     G4VPhysicalVolume* fCagewall5_phys;
     G4VPhysicalVolume* fCagewall6_phys;
     G4VPhysicalVolume* fMinicagewall_phys;
     G4VPhysicalVolume* fMinicagewall2_phys;
     G4VPhysicalVolume* fMinicagewall3_phys;
     G4VPhysicalVolume* fExtcagewall_phys;
     G4VPhysicalVolume* fExtcagewall2_phys;
     G4VPhysicalVolume* fExtboard_phys;
     G4VPhysicalVolume* fExtboard2_phys;
     G4VPhysicalVolume* fRing_phys;
     G4VPhysicalVolume* fFirstShield_botplate_phys;
     G4VPhysicalVolume* fFirstShield_side_phys;
     G4VPhysicalVolume* fFirstShield_topplate_phys;
     G4VPhysicalVolume* fFirstShield_topcyl_phys;
     G4VPhysicalVolume* fSecondShieldbotPlate_phys;
     G4VPhysicalVolume* fSecondShieldside_phys;
     G4VPhysicalVolume* fSecondShieldtopPlate_phys;
     G4VPhysicalVolume* fSecondShieldshield_topcyl_phys;
     G4VPhysicalVolume* fThirdShieldbotPlate_phys;
     G4VPhysicalVolume* fThirdShieldbotCone_phys;
     G4VPhysicalVolume* fThirdShieldtopCone_phys;
     G4VPhysicalVolume* fThirdShieldside_phys;
     G4VPhysicalVolume* fCryostatbotPlate_phys;
     G4VPhysicalVolume* fCryostatbotCone_phys;
     G4VPhysicalVolume* fCryostatbotSide_phys;
     G4VPhysicalVolume* fCryostatmidSide_phys;
     G4VPhysicalVolume* fCryostattopSide_phys;
     G4VPhysicalVolume* fCryostattopCone_phys;
     G4VPhysicalVolume* fAC_phys;
     G4VPhysicalVolume* fACIC_phys;
     G4VPhysicalVolume* fMesh_phys;
     G4VPhysicalVolume* fSphere_phys;
     G4VPhysicalVolume* fFiltercarrier_phys;
     G4VPhysicalVolume* fBSC_phys;
     G4VPhysicalVolume* fAlFilter_phys;
     G4VPhysicalVolume* fTrapezoid_A_phys;
     G4VPhysicalVolume* fTrapezoid_B_phys;
     G4VPhysicalVolume* fTrapezoid_C_phys;
     G4VPhysicalVolume* fTrapezoid_D_phys;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
