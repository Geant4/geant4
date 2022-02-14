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
#ifndef HGCALTBMATERIALS_HH
#define HGCALTBMATERIALS_HH

#include "G4Types.hh"
#include "G4String.hh"
#include "CLHEP/Units/SystemOfUnits.h"

#include <map>

class G4Box;
class G4LogicalVolume;
class G4Material;
class G4SubtractionSolid;

/// Construction of a hexagon solid
/// @param[in] name Name of the solid
/// @param[in] cellThickness Thickness of the hexagon
/// @param[in] cellSideLength Length of a side of the haxagon
/// @return Solid
G4SubtractionSolid *HexagonSolid(G4String name, G4double cellThickness,
                                 G4double cellSideLength);

/// Construction of a logical volume
/// @param[in] name Name of the solid
/// @param[in] cellThickness Thickness of the hexagon
/// @param[in] cellSideLength Length of a side of the haxagon
/// @param[in] material Material
/// @return Logical volume
G4LogicalVolume *HexagonLogical(G4String name, G4double cellThickness,
                                G4double cellSideLength, G4Material *material);

/**
 * @brief HGCal material and elements definitions
 *
 * Defines materials used for HGCal test-beam.
 * Creates logical volumes for every element that may be constructed and placed
 * in the test beam.
 *
 */

class HGCalTBMaterials {
public:
  /// Create HGCal materials
  HGCalTBMaterials();
  /// Set visualisation attributes
  void SetEventDisplayColorScheme();
  /// Get length of the beam line (world)
  inline G4double GetBeamLineLength() const { return fBeamLineLength; }
  /// Get transverse size of the beam line (world)
  inline G4double GetBeamLineXY() const { return fBeamLineXY; }
  /// Get pointer to the air material
  G4Material *GetAir() { return fMatAIR; }
  /// Place logical volume
  /// @param[in] aName Name of the logical volume
  /// @param[in,out] aZ0 Position in front of the logical volume, incremented by
  /// half the thickness for placement, and by another half thickness to return
  /// the position just behind the placed volume
  /// @param[in] aLogicMother Pointer to mother volume for placement
  void PlaceItemInLogicalVolume(std::string aName, G4double &aZ0,
                                G4LogicalVolume *aLogicMother);
  /// Get logical volume of silicon pixel (cell)
  G4LogicalVolume *GetSiPixelLogical() { return this->fSiPixelLogical; }
  /// Get logical volume of SiPM
  G4LogicalVolume *GetAHCALSiPMlogical() { return this->fAHCALSiPMlogical; }
  /// Get any logical volume by name
  inline const G4LogicalVolume *GetLogicalVolume(G4String aName) {
    return fLogicalVolumeMap[aName];
  };
  /// Get thickness of logical volume by name
  inline G4double GetThickness(std::string aName) {
    return fThicknessMap[aName];
  };

private:
  /// Define materisals used in HGCal test-beam
  void DefineMaterials();
  /// Define silicon wafer logical volume from silicon cells/pixels
  void DefineSiWaferAndCells();
  /// Define logical volumes for HGCal baseplates (CuW, Cu, PCB, Kapton layers)
  void DefineHGCalBaseplates();
  /// Define logical volumes for HGCal cases (Al, Steel)
  void DefineHGCalCases();
  /// Define logical volumes for HGCal electromagnetic part absorbers (Pb, Cu,
  /// W)
  void DefineHGCalEEAbsorbers();
  /// Define logical volumes for HGCal hadronic part absorbers (Cu, Fe)
  void DefineHGCalFHAbsorbers();
  /// Define logical volumes for AHCAL SiPM
  void DefineAHCALSiPM();
  /// Define logical volumes for AHCAL absorbers
  void DefineAHCALAbsorbers();
  /// Define logical volumes for beamline elements (MCP, scintillators, DWC)
  void DefineBeamLineElements();
  /// Length of the beam line
  G4double fBeamLineLength = 90 * CLHEP::m;
  /// Transverse dimension of the beam line
  G4double fBeamLineXY = 4 * CLHEP::m;
  /// Rotation angle of silicon hexagon
  G4double fAlpha;
  /// Side length of silicon cell hexagon
  G4double fSiPixelSideLength;
  /// Thickness of silicon wafer hexagon
  G4double fSiWaferThickness;
  /// Side length of silicon wafer haxagon
  G4double fSiWaferSideLength;
  /// Transverse size of AHCAL SiPM
  G4double fAHCALSiPMxy;
  /// Box representing AHCAL SiPM
  G4Box *fAHCALSiPMsolid;
  /// Map of volume name to its thickness
  std::map<G4String, G4double> fThicknessMap;
  /// Map of volume name to its logical volume
  std::map<G4String, G4LogicalVolume *> fLogicalVolumeMap;
  /// Map of volume name to counter of placed copies
  std::map<G4String, int> fCopyCounterMap;
  /// Materials
  G4Material *fMatVacuum;
  G4Material *fMatAIR;
  G4Material *fMatAr;
  G4Material *fMatAl;
  G4Material *fMatFe;
  G4Material *fMatGlass;
  G4Material *fMatSteel;
  G4Material *fMatPb;
  G4Material *fMatCu;
  G4Material *fMatW;
  G4Material *fMatSi;
  G4Material *fMatKAPTON;
  G4Material *fMatAu;
  G4Material *fMatPCB;
  G4Material *fMatQuartz;
  G4Material *fMatPolystyrene;
  G4Material *fMatCuW;
  G4Material *fMatC;
  G4Material *fMatH;
  G4Material *fMatO;
  G4Material *fMatMn;
  G4Material *fMatCr;
  G4Material *fMatNi;
  G4Material *fMatPolyethylene;
  G4Material *fMatFreon;
  G4Material *fMatScintillator;
  G4Material *fMatArCO2;
  G4Material *fMatCl;
  G4Material *fMatF;

  /// Logical volumes
  G4LogicalVolume *fSiPixelLogical;
  G4LogicalVolume *fSiWaferLogical;
  G4LogicalVolume *fCuWbaseplateLogical;
  G4LogicalVolume *fCuWbaseplate550umLogical;
  G4LogicalVolume *fCuWbaseplate610umLogical;
  G4LogicalVolume *fCuWbaseplate710umLogical;
  G4LogicalVolume *fCuBaseplateLogical;
  G4LogicalVolume *fCuBaseplate25umLogical;
  G4LogicalVolume *fCuBaseplate175umLogical;
  G4LogicalVolume *fPCBbaseplateLogical;
  G4LogicalVolume *fPCBbaseplateThinLogical;
  G4LogicalVolume *fKaptonLayerLogical;
  G4LogicalVolume *fAlCaseLogical;
  G4LogicalVolume *fAlCaseThickLogical;
  G4LogicalVolume *fAlChipLogical;
  G4LogicalVolume *fSteelCaseLogical;
  G4LogicalVolume *fSteelCaseThickLogical;
  G4LogicalVolume *fPbAbsorberEElogical;
  G4LogicalVolume *fFeAbsorberEElogical;
  G4LogicalVolume *fCuAbsorberEElogical;
  G4LogicalVolume *fWabsorberEElogical;
  G4LogicalVolume *fW2mmAbsorberEEDESY2018Logical;
  G4LogicalVolume *fW4mmAbsorberEEDESY2018Logical;
  G4LogicalVolume *fCuAbsorberFHlogical;
  G4LogicalVolume *fFeAbsorberFHlogical;
  G4LogicalVolume *fAHCALSiPMlogical;
  G4LogicalVolume *fAHCALSiPM2x2HUBlogical;
  G4LogicalVolume *fAlAbsorberAHCALlogical;
  G4LogicalVolume *fPCBAHCALlogical;
  G4LogicalVolume *fFeAbsorberAHCALlogical;
  G4LogicalVolume *fScintillatorLogical;
  G4LogicalVolume *fScintillatorThinLogical;
  G4LogicalVolume *fMCPlogical;
  G4LogicalVolume *fDWClogical;
  G4LogicalVolume *fDWCgasLogical;
  G4LogicalVolume *fCK3logical;
};
#endif /* HGCALTBMATERIALS_HH */
