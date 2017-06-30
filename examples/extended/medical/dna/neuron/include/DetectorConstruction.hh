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
// $ID$
/// \file DetectorConstruction.hh 
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VSolid.hh"
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

const G4int kMT = 100000;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Region;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DetectorConstruction();
  virtual ~DetectorConstruction();
  virtual G4VPhysicalVolume* Construct();
  G4Region* GetTargetRegion() {return fpRegion;}    

// For neuron constructions!
  NeuronLoadDataFile * GetNeuronLoadParamz() {return fNeuronLoadParamz;};  
  const G4VPhysicalVolume* GetsomaPV(G4int i) {return fsomaPV[i];};
  const G4VPhysicalVolume* GetdendritePV(G4int i) {return fdendritePV[i];};
  const G4VPhysicalVolume* GetaxonPV(G4int i) {return faxonPV[i];};
  G4Material* GetTargetMaterial()  {return fpWaterMaterial;};
  
private:
  G4Material*        fpDefaultMaterial;
  G4Material*        fpWaterMaterial;
  G4Region*          fpRegion;

  void DefineMaterials(); 
  G4VPhysicalVolume* ConstructDetector();    
  
  // For neuron constructions!  
  NeuronLoadDataFile * fNeuronLoadParamz;    
  G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps    
  G4VSolid * fsomaS[kMT], *fdendriteS[kMT], *faxonS[kMT];      
  G4LogicalVolume * fsomaLV[kMT], *fdendriteLV[kMT], *faxonLV[kMT];
  G4VPhysicalVolume *  fsomaPV[kMT], *fdendritePV[kMT], *faxonPV[kMT];

  G4VisAttributes * fSomaColour;
  G4VisAttributes * fDendColour;
  G4VisAttributes * fAxonColour;
  G4VisAttributes * fSpineColour;    
  G4VisAttributes * fNeuronColour;  

public:

  G4int   GetnbSomacomp()   {return fnbSomacomp;} //
  G4double  GetMassSomacomp (G4int i) {return fMassSomacomp[i];} //
  G4double GetMassSomaTot () {return fMassSomaTot;} //
  G4ThreeVector GetPosSomacomp(G4int i) {return fPosSomacomp[i];} //
  
  G4int   GetnbDendritecomp()   {return fnbDendritecomp;} //
  G4double GetMassDendcomp (G4int i) {return fMassDendcomp[i];} //  
  G4double GetMassDendTot () {return fMassDendTot;}  //
  G4ThreeVector GetPosDendcomp(G4int i) {return fPosDendcomp[i];} //
  G4double  GetDistADendSoma (G4int i) {return fDistADendSoma[i];} //
  G4double  GetDistBDendSoma (G4int i) {return fDistBDendSoma[i];} //
  
  G4int   GetnbAxoncomp()   {return fnbAxoncomp;} //
  G4double  GetMassAxoncomp (G4int i) {return fMassAxoncomp[i];} //
  G4double GetMassAxonTot () {return fMassAxonTot;}  //
  G4ThreeVector GetPosAxoncomp(G4int i) {return fPosAxoncomp[i];} //
  G4double  GetDistAxonsoma (G4int i) {return fDistAxonsoma[i];} //

  G4double GetTotVolNeuron() {return fTotVolNeuron;}
  G4double GetTotSurfNeuron () {return fTotSurfNeuron;} 
  G4double GetTotMassNeuron () {return fTotMassNeuron;} //  
  G4double GetTotVolSlice() {return fTotVolSlice;} 
  G4double GetTotSurfSlice () {return fTotSurfSlice;} 
  G4double GetTotMassSlice() {return fTotMassSlice;}  //
  G4double GetTotVolMedium() {return fTotVolMedium;} 
  G4double GetTotSurfMedium() {return fTotSurfMedium;} //
  G4double GetTotMassMedium() {return fTotMassMedium;}  //
 
private:

  G4int fnbSomacomp; 
  G4int fnbDendritecomp;
  G4int fnbAxoncomp;
  
  G4ThreeVector * fPosSomacomp ; 
  G4double * fMassSomacomp ;
  G4double fMassSomaTot ;
  
  G4double * fDistADendSoma ;
  G4double * fDistBDendSoma ;
  G4double * fMassDendcomp ;
  G4double fMassDendTot ;
  G4ThreeVector * fPosDendcomp ; 

  G4double * fDistAxonsoma ;
  G4double * fMassAxoncomp ;
  G4double fMassAxonTot ;
  G4ThreeVector * fPosAxoncomp ;  
 
  G4double fTotVolNeuron ;
  G4double fTotSurfNeuron ;
  G4double fTotMassNeuron ;  
  G4double fTotVolSlice;
  G4double fTotSurfSlice ;
  G4double fTotMassSlice; 
  G4double fTotVolMedium ;
  G4double fTotSurfMedium ;
  G4double fTotMassMedium ; 

};
#endif
