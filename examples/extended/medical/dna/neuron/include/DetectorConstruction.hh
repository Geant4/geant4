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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// and papers
// M. Batmunkh et al. J Radiat Res Appl Sci 8 (2015) 498-507
// O. Belov et al. Physica Medica 32 (2016) 1510-1520 
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
//
/// \file DetectorConstruction.hh 
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
//
#include "G4PVParameterised.hh"
#include "NeuronLoadDataFile.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Region;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DetectorConstruction();
  ~DetectorConstruction() override;

  G4VPhysicalVolume* Construct() override;
  G4Region* GetTargetRegion() const { return fpRegion; }    

  NeuronLoadDataFile * GetNeuronLoadParamz() const {return fNeuronLoadParamz;};  
  const G4VPhysicalVolume* GetsomaPV(G4int i) const {return fsomaPV[i];};
  const G4VPhysicalVolume* GetdendritePV(G4int i) const {return fdendritePV[i];};
  const G4VPhysicalVolume* GetaxonPV(G4int i) const {return faxonPV[i];};
  G4Material* GetTargetMaterial() const {return fpWaterMaterial;};

  G4int GetnbSomacomp() const {return fnbSomacomp;} //
  G4double GetMassSomacomp(G4int i) const {return fMassSomacomp[i];} //
  G4double GetMassSomaTot() const {return fMassSomaTot;} //
  G4ThreeVector GetPosSomacomp(G4int i) const {return fPosSomacomp[i];} //
  
  G4int GetnbDendritecomp() const {return fnbDendritecomp;} //
  G4double GetMassDendcomp(G4int i) const {return fMassDendcomp[i];} //  
  G4double GetMassDendTot() const {return fMassDendTot;}  //
  G4ThreeVector GetPosDendcomp(G4int i) const {return fPosDendcomp[i];} //
  G4double  GetDistADendSoma(G4int i) const {return fDistADendSoma[i];} //
  G4double  GetDistBDendSoma(G4int i) const {return fDistBDendSoma[i];} //
  
  G4int GetnbAxoncomp() const {return fnbAxoncomp;} //
  G4double GetMassAxoncomp(G4int i) const {return fMassAxoncomp[i];} //
  G4double GetMassAxonTot() const {return fMassAxonTot;}  //
  G4ThreeVector GetPosAxoncomp(G4int i) const {return fPosAxoncomp[i];} //
  G4double GetDistAxonsoma(G4int i) const {return fDistAxonsoma[i];} //

  G4double GetTotVolNeuron() const {return fTotVolNeuron;}
  G4double GetTotSurfNeuron() const {return fTotSurfNeuron;} 
  G4double GetTotMassNeuron() const {return fTotMassNeuron;} //  
  G4double GetTotVolSlice() const {return fTotVolSlice;} 
  G4double GetTotSurfSlice() const {return fTotSurfSlice;} 
  G4double GetTotMassSlice() const {return fTotMassSlice;}  //
  G4double GetTotVolMedium() const {return fTotVolMedium;} 
  G4double GetTotSurfMedium() const {return fTotSurfMedium;} //
  G4double GetTotMassMedium() const {return fTotMassMedium;}  //
  
private:

  void DefineMaterials(); 

  G4Material* fpWaterMaterial{nullptr};
  G4Material* fpWorldMaterial{nullptr};
  G4Region* fpRegion{nullptr};
  
  // For neuron constructions!  
  NeuronLoadDataFile* fNeuronLoadParamz;    
  G4bool fCheckOverlaps{false};

  std::vector<G4VSolid*> fsomaS; 
  std::vector<G4VSolid*> fdendriteS;
  std::vector<G4VSolid*> faxonS;      
  std::vector<G4LogicalVolume*> fsomaLV;
  std::vector<G4LogicalVolume*> fdendriteLV;
  std::vector<G4LogicalVolume*> faxonLV;
  std::vector<G4VPhysicalVolume*> fsomaPV;
  std::vector<G4VPhysicalVolume*> fdendritePV;
  std::vector<G4VPhysicalVolume*> faxonPV;

  G4VisAttributes * fSomaColour{nullptr};
  G4VisAttributes * fDendColour{nullptr};
  G4VisAttributes * fAxonColour{nullptr};
  G4VisAttributes * fSpineColour{nullptr};
  G4VisAttributes * fNeuronColour{nullptr};

  G4int fnbSomacomp{0}; 
  G4int fnbDendritecomp{0};
  G4int fnbAxoncomp{0};
  
  std::vector<G4ThreeVector> fPosSomacomp; 
  std::vector<G4double> fMassSomacomp;
  
  std::vector<G4double> fDistADendSoma;
  std::vector<G4double> fDistBDendSoma;
  std::vector<G4double> fMassDendcomp;
  std::vector<G4ThreeVector> fPosDendcomp; 

  std::vector<G4double> fDistAxonsoma;
  std::vector<G4double> fMassAxoncomp;
  std::vector<G4ThreeVector> fPosAxoncomp;  

  G4double fMassSomaTot{0.0};
  G4double fMassDendTot{0.0};
  G4double fMassAxonTot{0.0};
 
  G4double fTotVolNeuron{0.0};
  G4double fTotSurfNeuron{0.0};
  G4double fTotMassNeuron{0.0};
  G4double fTotVolSlice{0.0};
  G4double fTotSurfSlice{0.0};
  G4double fTotMassSlice{0.0};
  G4double fTotVolMedium{0.0};
  G4double fTotSurfMedium{0.0};
  G4double fTotMassMedium{0.0};

};

#endif
