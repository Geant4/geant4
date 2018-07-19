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
// $Id: 
// 
/// \file NeuronLoadDataFile.hh
/// \brief Implementation of the NeuronLoadDataFile class

#ifndef NeuronLoadDataFile_H
#define NeuronLoadDataFile_H 1

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
  virtual   
   ~NeuronLoadDataFile();  

  void SingleNeuronSWCfile(const G4String& filename);
  void NeuralNetworkDATAfile(const G4String& filename); 
  
// position, rotation of solids   
    void ComputeTransformation (const G4int copyNo,
                                G4VPhysicalVolume* physVol) const;
// Solid options: sphere or cylinder ...
   //G4VSolid* ComputeSolid (const G4int copyNo, 
 //   G4VPhysicalVolume* physiVol);
// and ... for solids!
    void ComputeDimensions(G4Tubs& cylinderComp, const G4int copyNo,
                                   const G4VPhysicalVolume*) const ;
    void ComputeDimensions(G4Sphere& , const G4int ,
                                   const G4VPhysicalVolume*) const {}         
    void ComputeDimensions(G4Ellipsoid& , const G4int ,
                                   const G4VPhysicalVolume*) const {}     
    void ComputeDimensions(G4Box&, const G4int, 
                                   const G4VPhysicalVolume*) const {}      
    void ComputeDimensions(G4Cons&,const G4int, 
                                   const G4VPhysicalVolume*) const {}
    void ComputeDimensions(G4Hype &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {} 
    void ComputeDimensions(G4Trd &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {} 
    void ComputeDimensions(G4Trap &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {} 
    void ComputeDimensions(G4Orb &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Torus &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Para &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Polycone &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
    void ComputeDimensions(G4Polyhedra &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
           
  G4double GetwidthB () {return fwidthB;} 
  G4double GetheightB () {return fheightB;} 
  G4double GetdepthB () {return fdepthB;} 
  G4double GetdiagnlLength () {return fdiagnlLength;} 
  G4double GetshiftX () {return fshiftX;} 
  G4double GetshiftY () {return fshiftY;} 
  G4double GetshiftZ () {return fshiftZ;} 
  G4double GetTypeN (G4int i) {return fTypeN[i];}
  
  G4int   GetnbSomacomp()   {return fnbSomacomp;} 
  G4double  GetMassSomacomp (G4int i) {return fMassSomacomp[i];}
  G4double GetMassSomaTot () {return fMassSomaTot;} 
  G4ThreeVector GetPosSomacomp(G4int i) {return fPosSomacomp[i];}
  G4double GetRadSomacomp(G4int i) {return fRadSomacomp[i];}
  
  G4int   GetnbDendritecomp()   {return fnbDendritecomp;} 
  G4double GetMassDendcomp (G4int i) {return fMassDendcomp[i];}  
  G4double GetMassDendTot () {return fMassDendTot;} 
  G4ThreeVector GetPosDendcomp(G4int i) {return fPosDendcomp[i];}
  G4double GetRadDendcomp(G4int i) {return fRadDendcomp[i];}
  G4double GetHeightDendcomp(G4int i) {return fHeightDendcomp[i];}
  G4double  GetDistADendSoma (G4int i) {return fDistADendSoma[i];}
  G4double  GetDistBDendSoma (G4int i) {return fDistBDendSoma[i];}
  G4RotationMatrix GetRotDendcomp(G4int i) {return fRotDendcomp[i];}
  
  G4int   GetnbAxoncomp()   {return fnbAxoncomp;} 
  G4double  GetMassAxoncomp (G4int i) {return fMassAxoncomp[i];}
  G4double GetMassAxonTot () {return fMassAxonTot;} 
  G4ThreeVector GetPosAxoncomp(G4int i) {return fPosAxoncomp[i];}
  G4double GetRadAxoncomp(G4int i) {return fRadAxoncomp[i];}
  G4double GetHeightAxoncomp(G4int i) {return fHeightAxoncomp[i];}
  G4double  GetDistAxonsoma (G4int i) {return fDistAxonsoma[i];}
  G4RotationMatrix GetRotAxoncomp(G4int i) {return fRotAxoncomp[i];}
 
  G4int   GetnbSpinecomp()   {return fnbSpinecomp;} 
  G4double  GetMassSpinecomp (G4int i) {return fMassSpinecomp[i];}
  G4double GetMassSpineTot () {return fMassSpineTot;} 
  G4ThreeVector GetPosSpinecomp(G4int i) {return fPosSpinecomp[i];}
  G4double GetRadSpinecomp(G4int i) {return fRadSpinecomp[i];}
  G4double GetHeightSpinecomp(G4int i) {return fHeightSpinecomp[i];}
  G4double  GetDistSpinesoma (G4int i) {return fDistSpinesoma[i];}
  G4RotationMatrix GetRotSpinecomp(G4int i) {return fRotSpinecomp[i];}   
  
  G4int   GetnbNeuroncomp()   {return fnbNeuroncomp;} 

  G4double GetTotVolNeuron() {return fTotVolNeuron;}
  G4double GetTotSurfNeuron () {return fTotSurfNeuron;} 
  G4double GetTotMassNeuron () {return fTotMassNeuron;}   
  G4double GetTotVolSlice() {return fTotVolSlice;} 
  G4double GetTotSurfSlice () {return fTotSurfSlice;} 
  G4double GetTotMassSlice() {return fTotMassSlice;}  
  G4double GetTotVolMedium() {return fTotVolMedium;} 
  G4double GetTotSurfMedium() {return fTotSurfMedium;} 
  G4double GetTotMassMedium() {return fTotMassMedium;}  
 
  G4VisAttributes GetSomaColour() {return fSomaColour;}
  G4VisAttributes GetDendColour() {return fDendColour;}
  G4VisAttributes GetAxonColour() {return fAxonColour;}
  G4VisAttributes GetSpineColour() {return fSpineColour;}
  G4VisAttributes GetNeuronColour() {return fNeuronColour;} 

  private:

  //! NEURON filename
  G4String fNeuronFileNameSWC;  
  G4String fNeuronFileNameDATA;  

  G4int fnbSomacomp; 
  G4int fnbDendritecomp;
  G4int fnbAxoncomp;
  G4int fnbSpinecomp;
  G4int fnbNeuroncomp;
  
  G4int * fnNn ;
  G4int * fpNn ;
  //G4int * ftypeC;
  G4int * fnNd ;
  G4int * fpNd ;
  G4int * fnNa ;
  G4int * fpNa ;
  G4int * fTypeN ;
  
  G4double fshiftX, fshiftY, fshiftZ; // shift in oder to center VOLUME!
  G4double fwidthB;
  G4double fheightB;
  G4double fdepthB;
  G4double fdiagnlLength;  // diagonal and diameter
  
  G4ThreeVector * fPosSomacomp ; 
  G4double * fRadSomacomp ;
  G4double * fMassSomacomp ;
  G4double fMassSomaTot ;
  
  G4double * fRadDendcomp ;
  G4double * fDistADendSoma ;
  G4double * fDistBDendSoma ;
  G4double * fHeightDendcomp ;
  G4double * fMassDendcomp ;
  G4double fMassDendTot ;
  G4ThreeVector * fPosDendcomp ; // VOXEL COORDINATES OF DENDRITES
  G4RotationMatrix * fRotDendcomp ; // RotationMatrix with Inverse

  G4double * fRadAxoncomp ;
  G4double * fHeightAxoncomp ;
  G4double * fDistAxonsoma ;
  G4double * fMassAxoncomp ;
  G4double fMassAxonTot ;
  G4ThreeVector * fPosAxoncomp ;  // VOXEL COORDINATES OF AXON
  G4RotationMatrix *  fRotAxoncomp ;

  G4double * fRadSpinecomp ;
  G4double * fHeightSpinecomp ;
  G4double * fDistSpinesoma ;
  G4double * fMassSpinecomp ;
  G4double fMassSpineTot ; 
  G4ThreeVector * fPosSpinecomp ;  // VOXEL COORDINATES OF SPINE
  G4RotationMatrix *  fRotSpinecomp ;
 
  G4double * fRadNeuroncomp ;
  G4double * fHeightNeuroncomp ;
  G4double * fDistNeuronsoma ;
  G4double * fMassNeuroncomp ;
  //G4double fMassNeuronTot ;
  G4ThreeVector * fPosNeuroncomp ;  // VOXEL COORDINATES OF Neuron
  G4RotationMatrix *  fRotNeuroncomp ;
  
  G4double fTotVolNeuron ;
  G4double fTotSurfNeuron ;
  G4double fTotMassNeuron ;  
  G4double fTotVolSlice;
  G4double fTotSurfSlice ;
  G4double fTotMassSlice; 
  G4double fTotVolMedium ;
  G4double fTotSurfMedium ;
  G4double fTotMassMedium ; 
 
  G4VisAttributes * fSomaColour;
  G4VisAttributes * fDendColour;
  G4VisAttributes * fAxonColour;
  G4VisAttributes * fSpineColour;    
  G4VisAttributes * fNeuronColour;  
  
  //G4Sphere* fsphereComp;
  //G4Tubs* fcylinderConp;

  //NeuronLoadMessenger * fpNeuronMessenger;

};

#endif


