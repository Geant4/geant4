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
// 
/// \file NeuronLoadDataFile.hh
/// \brief Implementation of the NeuronLoadDataFile class

#ifndef NeuronLoadDataFile_H
#define NeuronLoadDataFile_H 1

#include <vector>
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"

// Geant4 Constructive Solid Geometry (CSG)
#include "G4VPVParameterisation.hh"
#include "G4Box.hh"   // bounding volume
#include "G4Tubs.hh"   // axon, dendrite compartments
#include "G4Sphere.hh"   // soma compartments
#include "G4Ellipsoid.hh"  // soma compartments
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"   // stubby spine, filopodia
#include "G4EllipticalTube.hh"
#include "G4EllipticalCone.hh"
#include "G4Orb.hh"   // mushroom spine   
#include "G4Torus.hh"
#include "G4Para.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4Hype.hh"   // thin spine neck  
#include "G4Tet.hh"
#include "G4TwistedBox.hh"
#include "G4TwistedTrap.hh"
#include "G4TwistedTrd.hh"
#include "G4TwistedTubs.hh"

//class NeuronLoadMessenger;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class NeuronLoadDataFile 

// default - G4PVPlacement volume
// if use G4PVParameterised volume, 
// please enable following G4VPVParameterisation class!
//: public G4VPVParameterisation 

{ 
public:  
  NeuronLoadDataFile();
  ~NeuronLoadDataFile() = default;  

  void SingleNeuronSWCfile(const G4String& filename);
  void NeuralNetworkDATAfile(const G4String& filename); 
  
  // position, rotation of solids   
  void ComputeTransformation (const G4int copyNo,
			      G4VPhysicalVolume* physVol) const;

  void ComputeDimensions(G4Tubs& cylinderComp, const G4int copyNo,
			 const G4VPhysicalVolume*) const;
  G4double GetwidthB() const { return fwidthB; } 
  G4double GetheightB() const { return fheightB; } 
  G4double GetdepthB() const { return fdepthB; } 
  G4double GetdiagnlLength() const { return fdiagnlLength;} 
  G4double GetshiftX() const { return fshiftX; } 
  G4double GetshiftY() const { return fshiftY; } 
  G4double GetshiftZ() const { return fshiftZ; }
  G4double GetTypeN(G4int i) { return fTypeN[i]; }
  
  G4int GetnbSomacomp() const { return fnbSomacomp; } 
  G4double GetMassSomacomp(G4int i) const { return fMassSomacomp[i]; }
  G4double GetMassSomaTot() { return fMassSomaTot; } 
  G4ThreeVector GetPosSomacomp(G4int i) const { return fPosSomacomp[i]; }
  G4double GetRadSomacomp(G4int i) const { return fRadSomacomp[i]; }
  
  G4int GetnbDendritecomp() const { return fnbDendritecomp; } 
  G4double GetMassDendcomp(G4int i) const { return fMassDendcomp[i]; }  
  G4double GetMassDendTot() { return fMassDendTot; } 
  G4ThreeVector GetPosDendcomp(G4int i) const { return fPosDendcomp[i]; }
  G4double GetRadDendcomp(G4int i) const { return fRadDendcomp[i]; }
  G4double GetHeightDendcomp(G4int i) const { return fHeightDendcomp[i]; }
  G4double GetDistADendSoma(G4int i) const { return fDistADendSoma[i]; }
  G4double GetDistBDendSoma(G4int i) const { return fDistBDendSoma[i]; }
  G4RotationMatrix GetRotDendcomp(G4int i) const { return fRotDendcomp[i]; }
  
  G4int GetnbAxoncomp() const { return fnbAxoncomp; } 
  G4double GetMassAxoncomp (G4int i) const { return fMassAxoncomp[i]; }
  G4double GetMassAxonTot() const { return fMassAxonTot; } 
  G4ThreeVector GetPosAxoncomp(G4int i) const { return fPosAxoncomp[i]; }
  G4double GetRadAxoncomp(G4int i) const { return fRadAxoncomp[i]; }
  G4double GetHeightAxoncomp(G4int i) const { return fHeightAxoncomp[i]; }
  G4double GetDistAxonsoma(G4int i) const { return fDistAxonsoma[i]; }
  G4RotationMatrix GetRotAxoncomp(G4int i) const { return fRotAxoncomp[i]; }
 
  G4int GetnbSpinecomp() const { return fnbSpinecomp; } 
  G4double GetMassSpinecomp (G4int i) const { return fMassSpinecomp[i]; }
  G4double GetMassSpineTot() const { return fMassSpineTot; } 
  G4ThreeVector GetPosSpinecomp(G4int i) const { return fPosSpinecomp[i]; }
  G4double GetRadSpinecomp(G4int i) const { return fRadSpinecomp[i]; }
  G4double GetHeightSpinecomp(G4int i) const { return fHeightSpinecomp[i]; }
  G4double GetDistSpinesoma (G4int i) {return fDistSpinesoma[i];}
  G4RotationMatrix GetRotSpinecomp(G4int i) const { return fRotSpinecomp[i]; }  
  
  G4int GetnbNeuroncomp() const { return fnbNeuroncomp; }

  G4double GetTotVolNeuron() const { return fTotVolNeuron; }
  G4double GetTotSurfNeuron() const { return fTotSurfNeuron; } 
  G4double GetTotMassNeuron() const { return fTotMassNeuron; } 
  G4double GetTotVolSlice() const { return fTotVolSlice; } 
  G4double GetTotSurfSlice() const { return fTotSurfSlice; } 
  G4double GetTotMassSlice() const { return fTotMassSlice; }  
  G4double GetTotVolMedium() const { return fTotVolMedium; } 
  G4double GetTotSurfMedium() const { return fTotSurfMedium; } 
  G4double GetTotMassMedium() const { return fTotMassMedium; }  

private:

  //! NEURON filename
  G4String fNeuronFileNameSWC;  
  G4String fNeuronFileNameDATA;  

  G4int fnbSomacomp{0}; 
  G4int fnbDendritecomp{0};
  G4int fnbAxoncomp{0};
  G4int fnbSpinecomp{0};
  G4int fnbNeuroncomp{0};
  
  std::vector<G4int> fTypeN;
  
  // shift in oder to center VOLUME!
  G4double fshiftX{0.0};
  G4double fshiftY{0.0};
  G4double fshiftZ{0.0};
  G4double fwidthB{0.0};
  G4double fheightB{0.0};
  G4double fdepthB{0.0};
  G4double fdiagnlLength{0.0};  // diagonal and diameter

  G4double fMassSomaTot{0.0};
  G4double fMassDendTot{0.0};
  G4double fMassAxonTot{0.0};
  G4double fMassSpineTot{0.0};

  G4double fTotVolNeuron{0.0};
  G4double fTotSurfNeuron{0.0};
  G4double fTotMassNeuron{0.0};
  G4double fTotVolSlice{0.0};
  G4double fTotSurfSlice{0.0};
  G4double fTotMassSlice{0.0};
  G4double fTotVolMedium{0.0};
  G4double fTotSurfMedium{0.0};
  G4double fTotMassMedium{0.0};
  
  std::vector<G4ThreeVector> fPosSomacomp; 
  std::vector<G4double> fRadSomacomp;
  std::vector<G4double> fMassSomacomp;
  
  std::vector<G4double> fRadDendcomp;
  std::vector<G4double> fDistADendSoma;
  std::vector<G4double> fDistBDendSoma;
  std::vector<G4double> fHeightDendcomp;
  std::vector<G4double> fMassDendcomp;
  std::vector<G4ThreeVector> fPosDendcomp; // VOXEL COORDINATES OF DENDRITES
  std::vector<G4RotationMatrix> fRotDendcomp; // RotationMatrix with Inverse

  std::vector<G4double> fRadAxoncomp;
  std::vector<G4double> fHeightAxoncomp;
  std::vector<G4double> fDistAxonsoma;
  std::vector<G4double> fMassAxoncomp;
  std::vector<G4ThreeVector> fPosAxoncomp;  // VOXEL COORDINATES OF AXON
  std::vector<G4RotationMatrix> fRotAxoncomp;

  std::vector<G4double> fRadSpinecomp;
  std::vector<G4double> fHeightSpinecomp;
  std::vector<G4double> fDistSpinesoma;
  std::vector<G4double> fMassSpinecomp;
  std::vector<G4ThreeVector> fPosSpinecomp;  // VOXEL COORDINATES OF SPINE
  std::vector<G4RotationMatrix> fRotSpinecomp;
 
  std::vector<G4double> fRadNeuroncomp;
  std::vector<G4double> fHeightNeuroncomp;
  std::vector<G4double> fDistNeuronsoma;
  std::vector<G4double> fMassNeuroncomp;
  std::vector<G4ThreeVector> fPosNeuroncomp;  // VOXEL COORDINATES OF Neuron
  std::vector<G4RotationMatrix> fRotNeuroncomp;
};

#endif


