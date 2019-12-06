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
//  Author: F. Poignant, floriane.poignant@gmail.com
//
/// file STCyclotronDetectorConstruction.cc

#ifndef STCyclotronDetectorConstruction_h
#define STCyclotronDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include <fstream>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Region;
class G4Tubs;
class G4Material;
class STCyclotronDetectorMessenger;
class G4Element;

/// Detector construction class to define materials and geometry.

class STCyclotronDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  STCyclotronDetectorConstruction();
  ~STCyclotronDetectorConstruction();
    
  G4VPhysicalVolume* Construct();
  void ConstructSDandField();
  
  void SetTargetDiameter(G4double );
  void SetTargetIsotopeName(G4String );
  void SetTargetIsotopeZ(G4double );
  void SetTargetIsotopeN(G4int );
  void SetTargetIsotopeA(G4double );
  void SetTargetElementName(G4String );
  void SetTargetElementSymbole(G4String );
  void SetTargetElementNComponents(G4int );
  void SetTargetElementAbundance(G4double );
  void SetTargetMaterialDensity(G4double );
  void SetTargetMaterialNComponents(G4int );
  void SetTargetMaterialFractionMass(G4double );
  void SetTargetNaturalElement(G4String );
  void SetTargetNaturalMaterialFractionMass(G4double );
  G4bool UpdateMaterial();
  void SetTargetMaterial(G4String );

  void SetFoilIsotopeName(G4String );
  void SetFoilIsotopeZ(G4double );
  void SetFoilIsotopeN(G4int );
  void SetFoilIsotopeA(G4double );
  void SetFoilElementName(G4String );
  void SetFoilElementSymbole(G4String );
  void SetFoilElementNComponents(G4int );
  void SetFoilElementAbundance(G4double );
  void SetFoilMaterialDensity(G4double );
  void SetFoilMaterialNComponents(G4int );
  void SetFoilMaterialFractionMass(G4double );
  void SetFoilNaturalElement(G4String );
  void SetFoilNaturalMaterialFractionMass(G4double );
  G4bool UpdateFoilMaterial();
  void SetFoilMaterial(G4String );

  void SetTargetThickness(G4double );
  void SetFoilThickness(G4double );
  
  //Get methods
  inline G4double GetTargetPosition1(){return fLayer1_z_position_PART4 + 0.5*11.5 - fTarget_thickness;}
  inline G4double GetTargetPosition2(){return fLayer1_z_position_PART4 + 0.5*11.5;}
  inline G4double GetVolumeTarget(){return pi*fTarget_diameter*fTarget_diameter/4*fTarget_thickness;}
  inline G4double GetFoilPosition1(){return fZ_foil_position - 0.5*fFoil_thickness;}
  inline G4double GetTargetVolume(){return fTargetVolume;}
  inline G4double GetFoilVolume(){return fFoilVolume;}
  inline G4double GetFoilThickness(){return fFoil_thickness;}
  inline G4double GetTargetThickness(){return fTarget_thickness;}
  inline G4double GetTargetDiameter(){return fTarget_diameter;}
  
private:

  STCyclotronDetectorMessenger* fDetectorMessenger;

  //Messenger parameters
  G4double fTarget_diameter;
  std::vector<G4String> fIsotopeName;
  std::vector<G4double> fIsotopeZ;
  std::vector<G4int>    fIsotopeN;
  std::vector<G4double> fIsotopeA;
  std::vector<G4String> fElementName;
  std::vector<G4String> fElementSymbole;
  std::vector<G4int>    fElementNComponents;
  std::vector<G4double> fElementAbundance;
  std::vector<G4String> fNaturalElementName;
  std::vector<G4double> fNaturalMaterialFractionMass;
  G4double fDensity_target;
  G4int    fTarget_NComponents;
  std::vector<G4double> fMaterialFractionMass;

  std::vector<G4String> fIsotopeNameFoil;
  std::vector<G4double> fIsotopeZFoil;
  std::vector<G4int>    fIsotopeNFoil;
  std::vector<G4double> fIsotopeAFoil;
  std::vector<G4String> fElementNameFoil;
  std::vector<G4String> fElementSymboleFoil;
  std::vector<G4int>    fElementNComponentsFoil;
  std::vector<G4double> fElementAbundanceFoil;
  std::vector<G4String> fNaturalElementNameFoil;
  std::vector<G4double> fNaturalMaterialFractionMassFoil;
  G4double fDensity_foil;
  G4int    fFoil_NComponents;
  std::vector<G4double> fMaterialFractionMassFoil;

  G4double fTarget_thickness;
  G4double fFoil_thickness;

  //Parameters that are used/modified in the set methods
  //When modifying the target parameters
 
  //Material
  G4Material* fTarget_Material;
  G4Material* fFoil_Material;
  //Foil
  G4double fZ_foil_position;
  G4Tubs* fSolidFoil;
  G4LogicalVolume* fLogicFoil;
  G4VPhysicalVolume* fPhysFoil;
  //WORLD
  G4LogicalVolume* fLogicWorld;
  //PART 3
  G4double fLayer_z_position_PART3;
  G4VPhysicalVolume* fPhysLayer_PART3;
  G4VPhysicalVolume* fPhysTube_PART3;
  //PART 4
  G4double fTube_outerRadius_PART4;
  G4double fTube_length_PART4;
  G4double fLayer_z_position_PART4;
  G4VPhysicalVolume* fPhysTube_PART4;
  G4VPhysicalVolume* fPhysLayer_PART4;
  G4double fLayer1_z_position_PART4;
  G4VPhysicalVolume* fPhysLayer1_PART4;
  //Target
  G4LogicalVolume* fLogicTarget;
  G4double fTarget_z_position;
  G4Tubs* fSolidTarget;
  G4VPhysicalVolume* fPhysTarget;
  //PART 5 
  G4double fLayer1_z_position_PART5;
  G4VPhysicalVolume* fPhysLayer1_PART5;
  G4double fLayer2_z_position_PART5;
  G4VPhysicalVolume* fPhysLayer2_PART5;
  G4double fLayer3_z_position_PART5;
  G4VPhysicalVolume* fPhysLayer3_PART5;
  
  G4Region* fRegionTarget;
  G4Region* fRegionFoil;

  G4double fTargetVolume;
  G4double fFoilVolume;
  
  std::ofstream fParametersSummary;

};
#endif
