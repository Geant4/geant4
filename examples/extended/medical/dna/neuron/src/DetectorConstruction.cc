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
/// \file DetectorConstruction.cc 
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4PhysicalConstants.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"
#include <algorithm>  
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  fNeuronLoadParamz = new NeuronLoadDataFile();   
  DefineMaterials();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fNeuronLoadParamz;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  // Water is defined from NIST material database
  G4NistManager * man = G4NistManager::Instance();
  fpWaterMaterial = man->FindOrBuildMaterial("G4_WATER");
  fpWorldMaterial = man->FindOrBuildMaterial("G4_Galactic");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4cout <<" ---- Begin of Neuron Construction! ---- " 
	 <<"\n"<<" =========================================================="<<G4endl;

  // =============================================== 
  // WORLD VOLUME - filled by default material
  // ===============================================

  // Dimensions of world volume are calculated as overall dimensions of neuron! 
  
  G4double worldSizeX;
  worldSizeX = fNeuronLoadParamz->GetdiagnlLength()*um;

  if (worldSizeX <= 0.0) {
    worldSizeX = 10.*cm;
  }
 
  G4double worldSizeY = worldSizeX;
  G4double worldSizeZ = worldSizeX;
  G4cout << " Side length of word volume is calculated : " 
         << worldSizeX/um <<" um"<< G4endl;  
  G4VSolid* worldS = new G4Box("World", worldSizeX/2, worldSizeY/2, worldSizeZ/2);

  G4LogicalVolume* worldLV = new G4LogicalVolume(worldS, fpWorldMaterial, "World");

  // Visualization attributes
  G4VisAttributes* worldVisAtt =
      new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.1)); //Gray
  worldVisAtt->SetVisibility(true);
  worldLV->SetVisAttributes(worldVisAtt);
  
  G4VPhysicalVolume* worldPV = new G4PVPlacement(0, //no rotation
                                  G4ThreeVector(),  //at (0,0,0)
                                  "World",     //its name
                                  worldLV,     //its logical volume
                                  0,       //its mother  volume
                                  false,      //no boolean operation
                                  0,       //copy number
                                  true); // checking overlaps forced to YES

  // ===============================================
  // HOMOGENEOUS MEDIUM - LARGE SPHERE VOLUME
  // ===============================================

  // Radius of water sphere calculated as overall dimensions of neuron.
  fTotMassMedium = fNeuronLoadParamz->GetTotMassMedium() ;
  fTotSurfMedium = fNeuronLoadParamz->GetTotSurfMedium() ;
  G4double RadiusMedium = fNeuronLoadParamz->GetdiagnlLength()*um / 2.; 
  G4cout << " Radius of homogeneous medium is calculated : " 
         << RadiusMedium/um <<" um"<< G4endl; 
  G4VSolid* mediumS = new G4Orb("Medium", RadiusMedium);   
  
  G4LogicalVolume* mediumLV = new G4LogicalVolume(mediumS,  fpWaterMaterial, "Medium");    
        
  // Visualization attributes
  G4VisAttributes* mediumVisAtt =
      new G4VisAttributes(G4Colour(0.0,1.0,0.0,0.02)); //Green
  //mediumVisAtt->SetForceSolid(true);
  //mediumVisAtt->SetForceWireframe (true);
  mediumVisAtt->SetForceLineSegmentsPerCircle(180);
  mediumVisAtt->SetVisibility(true);
  mediumLV->SetVisAttributes(mediumVisAtt);
  
  G4VPhysicalVolume* mediumPV = new G4PVPlacement(0,  
                    G4ThreeVector(), 
                    "Medium",  
                    mediumLV,  
                    worldPV,  
                    false,        
                    0,             
                    fCheckOverlaps);

  // ===============================================
  // TARGET - BOUNDING SLICE including NEURON 
  // ===============================================

  // Dimensions of bounding slice volume defined as overall measure of neuron
 
  G4double TargetSizeX = fNeuronLoadParamz->GetwidthB()*um; 
  G4double TargetSizeY = fNeuronLoadParamz->GetheightB()*um; 
  G4double TargetSizeZ = fNeuronLoadParamz->GetdepthB()*um; 
  fTotMassSlice = fNeuronLoadParamz->GetTotMassSlice() ;  
  G4cout << " Overall dimensions (um) of neuron morphology : " << "\n"
         << '\t'<< " width = " <<TargetSizeX/um<< " height = " << 
         TargetSizeY/um
         << " depth = " <<TargetSizeZ/um<<G4endl;  
  
  G4cout << " Volume (um3), surface (um2) and mass (ug) of Bounding Slice are"
	 << " calculated : " << "\n"
	 << '\t'<<fNeuronLoadParamz->GetTotVolSlice()/std::pow(um,3)<<"; "<<'\t'
	 << fNeuronLoadParamz->GetTotSurfSlice()/(um*um)
	 << "; "<<'\t'<<fNeuronLoadParamz->GetTotMassSlice()*1e6/g<< "\n "<<G4endl;   
 
  G4Box* boundingS =
    new G4Box("BoundingSlice", TargetSizeX/2.,TargetSizeY/2.,TargetSizeZ/2.);   

  G4LogicalVolume* boundingLV =
    new G4LogicalVolume(boundingS, fpWaterMaterial, "BoundingSlice");    

  // Visualization attributes with opacity!
  G4VisAttributes* TargetVisAtt =
      new G4VisAttributes(G4Colour(1.0,1.0,1.0,0.1)); 
  TargetVisAtt->SetForceSolid(true);
  TargetVisAtt->SetVisibility(true);
  boundingLV->SetVisAttributes(TargetVisAtt); 
  new G4PVPlacement(0,      
                    G4ThreeVector(),   
                    "BoundingSlice",     
                    boundingLV,   
                    mediumPV,     
                    false,          
                    0,             
      fCheckOverlaps);     

  // ===============================================
  // NEURON MORPHOLOGY
  // ===============================================

  G4cout <<  " Volume (um3), surface (um2) and mass(ug) of Neuron "
        << "are calculated : "<< "\n"
        << '\t'<<fNeuronLoadParamz->GetTotVolNeuron()/std::pow(um,3)<<"; "<<'\t'
        << fNeuronLoadParamz->GetTotSurfNeuron()/(um*um)
        <<"; "<<'\t'<<fNeuronLoadParamz->GetTotMassNeuron()*1e6/g<<G4endl;  
  fTotMassNeuron = fNeuronLoadParamz->GetTotMassNeuron();  
  G4cout <<  " Total number of compartments into Neuron : "
         <<fNeuronLoadParamz->GetnbNeuroncomp()<<G4endl; 
  G4cout << " Shift values (um) for Neuron translation are calculated : " << "\n"
	 << '\t' << " shiftX = " <<fNeuronLoadParamz->GetshiftX()<< " shiftY = " 
	 << fNeuronLoadParamz->GetshiftY()
	 << " shiftZ = " <<fNeuronLoadParamz->GetshiftZ()<< G4endl;

  // Soma in Violet with opacity   // 0.85,0.44,0.84
  fSomaColour = new G4VisAttributes;
  fSomaColour->SetColour(G4Colour(G4Colour(22/255. , 200/255. , 30/255.))); 
  fSomaColour->SetForceSolid(true); // true
  fSomaColour->SetVisibility(true);
  
  // Dendrites in Dark-Blue  
  fDendColour = new G4VisAttributes;
  fDendColour->SetColour(G4Colour(G4Colour(0.0, 0.0, 0.5)));
  fDendColour->SetForceSolid(true);
  //fDendColour->SetVisibility(true);

  // Axon in Maroon  
  fAxonColour = new G4VisAttributes;
  fAxonColour->SetColour(G4Colour(G4Colour(0.5, 0.0, 0.0))); 
  fAxonColour->SetForceSolid(true);
  fAxonColour->SetVisibility(true);

  // Spines in Dark-Green   
  fSpineColour = new G4VisAttributes;
  fSpineColour->SetColour(G4Colour(G4Colour(0.0 , 100/255. , 0.0)));
  fSpineColour->SetForceSolid(true);
  fSpineColour->SetVisibility(true);    

  // Whole neuron in semitransparent navy blue   
  fNeuronColour = new G4VisAttributes;
  fNeuronColour->SetColour(G4Colour(G4Colour(0.0,0.4,0.8,0.9)));
  fNeuronColour->SetForceSolid(true);
  fNeuronColour->SetVisibility(true); 
 
  // Placement volumes: G4examples/extended/parameterisations/gflash 

  // =======================================================================
  // Structure-1: Soma

  // Create Target G4Region and add logical volume
  // Active Geant4-DNA processes in this region 
  fpRegion = new G4Region("Soma"); 
  G4ProductionCuts* cuts = new G4ProductionCuts();
  G4double defCut = 1*nanometer;
  cuts->SetProductionCut(defCut,"gamma");
  cuts->SetProductionCut(defCut,"e-");
  cuts->SetProductionCut(defCut,"e+");
  cuts->SetProductionCut(defCut,"proton");
  fpRegion->SetProductionCuts(cuts);
  G4ThreeVector shift(fNeuronLoadParamz->GetshiftX(),
		      fNeuronLoadParamz->GetshiftY(),
		      fNeuronLoadParamz->GetshiftZ());

  G4int n = fNeuronLoadParamz->GetnbSomacomp();
  if (n <= 0) {
    G4cout <<" ---- Soma not found! ---- "<< G4endl; 
  }
  else {
    G4cout <<" ---- Soma for construction: ---- "<< G4endl; 
    G4cout <<  " Total number of compartments into Soma : " << n << G4endl;
    fnbSomacomp = n;
    fMassSomaTot = fNeuronLoadParamz->GetMassSomaTot();
    fMassSomacomp.resize(n, 0.0);
    fPosSomacomp.resize(n);
    fsomaS.resize(n, nullptr);
    fsomaLV.resize(n, nullptr);
    fsomaPV.resize(n, nullptr);

    for (G4int i=0; i<n; ++i) {
      fsomaS[i] = new G4Orb("Soma", fNeuronLoadParamz->GetRadSomacomp(i)* um);
      // you can change parameters of Soma with a single compartment
      G4ThreeVector pos = (fNeuronLoadParamz->GetPosSomacomp(i) - shift)*um;
      fsomaLV[i] = new G4LogicalVolume(fsomaS[i], fpWaterMaterial, "Soma");
      fsomaLV[i]->SetVisAttributes(fSomaColour);
      fsomaPV[i] = new G4PVPlacement(0,  // no rotation
				     pos,
				     fsomaLV[i],
				     "Soma",
				     boundingLV,
				     false,
				     i);
      fMassSomacomp[i] = fNeuronLoadParamz->GetMassSomacomp(i) ;
      fPosSomacomp[i] = fNeuronLoadParamz->GetPosSomacomp(i) ;
      fpRegion->AddRootLogicalVolume(fsomaLV[i]); 
    }
  }  

  // ======================================================================= 
  // Structure-2: Dendrites

  n = fNeuronLoadParamz->GetnbDendritecomp();
  if (n <= 0) {
    G4cout <<" ---- Dendrites not found! ---- "<< G4endl; 
  }
  else {
    fnbDendritecomp = n;
    G4cout << " ---- Dendrites for construction: ---- " << G4endl; 
    G4cout << " Total number of compartments into Dendrites: "  << n <<G4endl; 

    // Active Geant4-DNA processes in this region 
    auto region = new G4Region("Dendrites");
    region->SetProductionCuts(cuts);

    fMassDendTot = fNeuronLoadParamz->GetMassDendTot();
    fMassDendcomp.resize(n, 0.0);
    fDistADendSoma.resize(n, 0.0);
    fDistBDendSoma.resize(n, 0.0);
    fPosDendcomp.resize(n);
    fdendriteS.resize(n, nullptr);
    fdendriteLV.resize(n, nullptr);
    fdendritePV.resize(n, nullptr);
    for (G4int i=1; i<n; ++i) {  
      fdendriteS[i] = new G4Tubs("Dendrites",
				 0., fNeuronLoadParamz->GetRadDendcomp(i)*um, 
				 fNeuronLoadParamz->GetHeightDendcomp(i)*um/2., 
				 0., 2.*pi ); 
      fdendriteLV[i] = new G4LogicalVolume(fdendriteS[i], fpWaterMaterial, "Dendrites");
      fdendriteLV[i]->SetVisAttributes(fDendColour);
 
      G4ThreeVector pos = (fNeuronLoadParamz->GetPosDendcomp(i) - shift)*um;
      // rotation checking with function ComputeTransformation!
      // RotationMatrix with Inverse
      fdendritePV[i] = new G4PVPlacement(G4Transform3D(fNeuronLoadParamz->GetRotDendcomp(i), pos),
					 fdendriteLV[i], "Dendrites",
					 boundingLV, false, i);
      fMassDendcomp[i] = fNeuronLoadParamz->GetMassDendcomp(i);
      fPosDendcomp[i] = fNeuronLoadParamz->GetPosDendcomp(i);
      fDistADendSoma[i] = fNeuronLoadParamz->GetDistADendSoma(i);
      fDistBDendSoma[i] = fNeuronLoadParamz->GetDistBDendSoma(i);
      region->AddRootLogicalVolume(fdendriteLV[i]);
    }
  }
 
  // =======================================================================
  // Structure-3: Axon 

  n = fNeuronLoadParamz->GetnbAxoncomp();
  if (n <= 0) {
    G4cout <<" ---- Axon not found! ---- "<< G4endl; 
  }
  else {
    G4cout <<" ---- Axon for construction: ---- "<< G4endl;   
    G4cout <<  " Total number of compartments into Axon : " << n << G4endl; 
    // Active Geant4-DNA processes in this region 
    auto region = new G4Region("Axon"); 
    region->SetProductionCuts(cuts);

    fnbAxoncomp = n;
    fMassAxonTot = fNeuronLoadParamz->GetMassAxonTot();
    fMassAxoncomp.resize(n, 0.0);
    fDistAxonsoma.resize(n, 0.0);
    fPosAxoncomp.resize(n);
    faxonS.resize(n, nullptr);
    faxonLV.resize(n, nullptr);
    faxonPV.resize(n, nullptr);

    for (G4int i=1; i<n; ++i) {
      faxonS[i] = new G4Tubs("Axon",
			     0., fNeuronLoadParamz->GetRadAxoncomp(i)*um, 
			     fNeuronLoadParamz->GetHeightAxoncomp(i)*um/2., 
			     0., 2.*pi ); 
      faxonLV[i] = new G4LogicalVolume(faxonS[i], fpWaterMaterial, "Axon");
      faxonLV[i]->SetVisAttributes(fAxonColour);
      // RotationMatrix with Inverse
      G4ThreeVector pos = (fNeuronLoadParamz->GetPosAxoncomp(i) - shift)*um;
      faxonPV[i] = new G4PVPlacement(G4Transform3D(fNeuronLoadParamz->GetRotAxoncomp(i), pos),
				     faxonLV[i], "Axon",
				     boundingLV, false, i);
      fMassAxoncomp[i] = fNeuronLoadParamz->GetMassAxoncomp(i);
      fPosAxoncomp[i] = fNeuronLoadParamz->GetPosAxoncomp(i);
      fDistAxonsoma[i] = fNeuronLoadParamz->GetDistAxonsoma(i);
      region->AddRootLogicalVolume(faxonLV[i]);
    }
  }   
  // =======================================================================
  // Structure-4: Spines 
  if (fNeuronLoadParamz->GetnbSpinecomp()==0)
  {
    G4cout <<" ---- Spines not found! ---- "<< G4endl; 
  }
  else {
    G4cout <<" ---- Spines for construction: ---- "<< G4endl; 
    G4cout <<  " Total number of compartments into Spines : "
	   << fNeuronLoadParamz->GetnbSpinecomp() << G4endl; 
  }
 
  G4cout <<"\n ---- End of Neuron Construction! ---- " 
      << "\n ========================================================== \n"
      << G4endl; 
  
  // =======================================================================
  // Active Geant4-DNA processes in BoundingSlice with whole neuron structure
  //fpRegion = new G4Region("BoundingSlice"); 
  //fpRegion->SetProductionCuts(cuts);
  // fpRegion->AddRootLogicalVolume(boundingLV); 

  return worldPV;
}
