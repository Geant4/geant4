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
/// \file NeuronLoadDataFile.cc
/// \brief Implementation of the NeuronLoadDataFile class

#include "NeuronLoadDataFile.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"
#include "G4ios.hh"
#include <algorithm>  
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include "globals.hh"
#include "CommandLineParser.hh"
#include "G4UImanager.hh"

using namespace std;
using namespace G4DNAPARSER;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
NeuronLoadDataFile::NeuronLoadDataFile()
{
  Command* commandLine(nullptr);
 
  // 1. Load single neuron morphology and obtain parameters.
  // Default SWC file name of neuron
  fNeuronFileNameSWC = "GranuleCell-Nr2.CNG.swc";  

  // 2. Load neural network and obtain parameters.
  // Default prepared data filename of neural network with single/multi-layer.
  // Small network of 10 pyramidal neurons with single layer
  fNeuronFileNameDATA = "NeuralNETWORK.dat"; 

  // Load/change SWC or DAT as "CommandLineParser" class
  if((commandLine=CommandLineParser::GetParser()->GetCommandIfActive("-swc"))) {
    fNeuronFileNameSWC = G4String(commandLine->GetOption());
  }
  if ((commandLine = CommandLineParser::GetParser()->GetCommandIfActive("-network"))) {
    auto ss = G4String(commandLine->GetOption());
    if ("" != ss) { fNeuronFileNameDATA = ss; }
    NeuralNetworkDATAfile(fNeuronFileNameDATA); 
  }
  else {
    SingleNeuronSWCfile(fNeuronFileNameSWC);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeuronLoadDataFile::SingleNeuronSWCfile(const G4String& filename)
{
  // -----------
  // 12 November 2012 - code created
  // -------------------------------------------------------------------
  // November 2012: First model of neuron[*] adapted into Geant4 microdosimetry  
  //                from Claiborne`s database[**] by M. Batmunkh.
  // February 2013: Loading SWC file from NeuronMorpho.Org[***] 
  //                suggested by L. Bayarchimeg.
  // [*] http://lt-jds.jinr.ru/record/62124/files/lrb_e_2012.pdf
  // [**] http://www.utsa.edu/claibornelab/
  // [***] http://neuromorpho.org
  // -------------------------------------------------------------------

  G4String sLine = "";
  std::ifstream infile;
  infile.open(filename.c_str());
  if (!infile.is_open()) {
    G4ExceptionDescription ed;
    ed << "Datafile " << filename << " is not opened!";
    G4Exception("NeuronLoadDataFile::SingleNeuronSWCfile()","dna014",
		FatalException, ed, "Check file path");
    return;
  } 
  G4int nrows,nlines;
  nrows=0; nlines=0;
  for(;;) {
    getline(infile, sLine);
    if (infile.eof()) { break; }
    ++nrows;
  }
  infile.close();
  
  G4cout << "NeuronLoadDataFile::SingleNeuronSWCfile: number of strings="<< nrows 
	 << " in file " << filename << G4endl;  
  
  infile.open(filename.c_str());
 
  G4double TotVolSoma, TotVolDend, TotVolAxon, TotVolSpine;
  TotVolSoma=TotVolDend=TotVolAxon=TotVolSpine=0.;
  G4double TotSurfSoma, TotSurfDend, TotSurfAxon, TotSurfSpine;
  TotSurfSoma=TotSurfDend=TotSurfAxon=TotSurfSpine=0.;
  G4int nNcomp;    // current index of neuronal compartment
  G4int typeNcomp; // type of neuron structures: soma, axon, dendrite, etc. 
  G4double x,y,z;  // cartesian coordinates of each compartment in micrometer
  G4double radius; // radius of each compartment in micrometer
  G4int pNcomp;    // linked compartment, indicates branch points of dendrites  
  G4double minX,minY,minZ;
  G4double maxX,maxY,maxZ; 
  G4double maxRad = -1e+09;
  minX=minY=minZ=1e+09;
  maxX=maxY=maxZ=-1e+09;
  G4double density = 1.0 * (g/cm3) ; // water medium 
  G4double Piconst = (4.0/3.0)*pi ; 

  fMassSomacomp.resize(nrows, 0);
  fPosSomacomp.resize(nrows);
  fRadSomacomp.resize(nrows, 0); 
  std::vector<G4ThreeVector> PosDendcomp(nrows);
  fRadDendcomp.resize(nrows, 0);
  fHeightDendcomp.resize(nrows, 0);
  fMassDendcomp.resize(nrows, 0);
  fDistADendSoma.resize(nrows, 0);
  fDistBDendSoma.resize(nrows, 0);
  fPosDendcomp.resize(nrows);
  fRotDendcomp.resize(nrows);
  std::vector<G4ThreeVector> PosAxoncomp(nrows);
  fRadAxoncomp.resize(nrows, 0);
  fHeightAxoncomp.resize(nrows, 0);
  fMassAxoncomp.resize(nrows, 0);
  fDistAxonsoma.resize(nrows, 0);
  fPosAxoncomp.resize(nrows);
  fRotAxoncomp.resize(nrows);
  fMassSpinecomp.resize(nrows, 0);
  fPosSpinecomp.resize(nrows);
  fRadSpinecomp.resize(nrows, 0); 
  fRadNeuroncomp.resize(nrows, 0);
  fHeightNeuroncomp.resize(nrows, 0);
  fDistNeuronsoma.resize(nrows, 0);
  fPosNeuroncomp.resize(nrows);
  fRotNeuroncomp.resize(nrows);
  fPosNeuroncomp.resize(nrows);
  fRadNeuroncomp.resize(nrows, 0); 
  fTypeN.resize(nrows, 0);
  G4ThreeVector base;
 
  // to read datafile containing numbers, alphabets and symbols..,
  for (;;) {
    getline(infile, sLine);
    if (infile.eof()) { break; }
    if ("#" == sLine.substr(0, 1)) { continue; };

    std::istringstream form(sLine);
    form >> nNcomp >> typeNcomp >> x >> y >> z >> radius >> pNcomp;
    /*
	G4cout << "NeuronLoadDataFile::SingleNeuronSWCfile: typeNcomp="
	       << typeNcomp << " nNcomp=" << nNcomp << " pNcomp=" << pNcomp << " N1="
	       << fnbSomacomp << " N2=" << fnbAxoncomp << " N3=" << fnbDendritecomp << G4endl;
    */
    // =======================================================================
    // to find the largest and the smallest values of compartment positions
    // for parameters of bounding slice, sphere medium and shift of neuron.
    if (minX > x) minX = x;
    if (minY > y) minY = y;
    if (minZ > z) minZ = z;
    if (maxX < x) maxX = x;
    if (maxY < y) maxY = y;
    if (maxZ < z) maxZ = z; 
    // max diameter of compartments 
    if (maxRad < radius) maxRad = radius; 
  
    // =======================================================================
    // Soma compartments represented as Sphere or Ellipsoid solid
    if (typeNcomp == 1) { 
      //  Sphere volume and surface area
      G4double VolSomacomp = Piconst*pow(radius*um, 3);
      TotVolSoma = TotVolSoma + VolSomacomp;
      G4double SurSomacomp = 3.*Piconst*pow(radius*um, 2);
      TotSurfSoma = TotSurfSoma + SurSomacomp;
      fMassSomacomp[fnbSomacomp] = density*VolSomacomp;
      fMassSomaTot = fMassSomaTot + fMassSomacomp[fnbSomacomp];  
      G4ThreeVector vSoma(x, y, z); 
      fPosSomacomp[fnbSomacomp] = vSoma; 
      fRadSomacomp[fnbSomacomp] = radius;
      if (0 == fnbSomacomp) {
	base = G4ThreeVector(fRadSomacomp[0], fRadSomacomp[0], fRadSomacomp[0]);
      }
      ++fnbSomacomp;
    }

    // =======================================================================
    // Apical and basal dendritic compartments represented as cylinderical solid
    if (typeNcomp == 3 || typeNcomp == 4) {
      G4ThreeVector vDend(x, y, z);
      // Position and Radius of compartments
      PosDendcomp[fnbDendritecomp] = vDend;
      fRadDendcomp[fnbDendritecomp] = radius;
      // To join two tracing points along the dendritic branches. 
      // To calculate length, center and rotation angles of each cylinder
      fPosDendcomp[fnbDendritecomp] = vDend; //translmDend;
      // delta of position A and position B of cylinder 
      G4ThreeVector dend;
      //primary dendritic branch should be connect with Soma
      if (0 == fnbDendritecomp) {
	dend = PosDendcomp[fnbDendritecomp] - fPosSomacomp[0] - base;
      }
      else {
	dend = PosDendcomp[fnbDendritecomp] - PosDendcomp[fnbDendritecomp - 1];
      }
      // Height of compartment
      G4double lengthDendcomp = dend.mag();
      fHeightDendcomp[fnbDendritecomp] = lengthDendcomp;

      // Distance from Soma
      G4ThreeVector dendDis = fPosSomacomp[0] - fPosDendcomp[fnbDendritecomp];
      if (typeNcomp == 3) fDistADendSoma[fnbDendritecomp] = dendDis.mag();
      if (typeNcomp == 4) fDistBDendSoma[fnbDendritecomp] = dendDis.mag();
   
      //  Cylinder volume and surface area
      G4double VolDendcomp = pi*pow(radius*um,2)*(lengthDendcomp*um);
      TotVolDend = TotVolDend + VolDendcomp;
      G4double SurDendcomp = 2.*pi*radius*um*(radius+lengthDendcomp)*um;
      TotSurfDend = TotSurfDend + SurDendcomp;
      fMassDendcomp[fnbDendritecomp] = density*VolDendcomp; 
      fMassDendTot = fMassDendTot + fMassDendcomp[fnbDendritecomp];   
   
      dend = dend.unit();
   
      // Euler angles of each compartment
      G4double theta_eulerDend = dend.theta();
      G4double phi_eulerDend = dend.phi();
      G4double psi_eulerDend = 0;

      //Rotation Matrix, Euler constructor build inverse matrix.
      G4RotationMatrix rotmDendInv  = G4RotationMatrix(phi_eulerDend+pi/2,
						       theta_eulerDend,
						       psi_eulerDend);
      fRotDendcomp[fnbDendritecomp] = rotmDendInv.inverse();
      ++fnbDendritecomp;
    }
  
    // =======================================================================
    // Axon compartments represented as cylinderical solid
    if (typeNcomp == 2 || typeNcomp == 7) {
      G4ThreeVector vAxon(x, y, z); 
      // Position and Radius of compartments
      PosAxoncomp[fnbAxoncomp] = vAxon;
      fRadAxoncomp[fnbAxoncomp] = radius;
      // To join two tracing points in loaded SWC data file. 
      // To calculate length, center and rotation angles of each cylinder   

      // delta of position A and position B of cylinder 
      G4ThreeVector Axon;
      //primary axon point should be connect with Soma
      if (0 == fnbAxoncomp) {
	Axon = PosAxoncomp[fnbAxoncomp] - fPosSomacomp[0] - base; 
      }
      else {
	Axon = PosAxoncomp[fnbAxoncomp] - PosAxoncomp[fnbAxoncomp - 1];
      }
      G4double lengthAxoncomp = Axon.mag();
      // Height of compartment
      fHeightAxoncomp[fnbAxoncomp] = lengthAxoncomp;
   
      // Distance from Soma
      G4ThreeVector AxonDis = fPosSomacomp[0] - fPosAxoncomp[fnbAxoncomp];
      fDistAxonsoma[fnbAxoncomp] = AxonDis.mag();
         
      //  Cylinder volume and surface area
      G4double VolAxoncomp = pi*pow(radius*um, 2)*(lengthAxoncomp*um);
      TotVolAxon = TotVolAxon + VolAxoncomp;
      G4double SurAxoncomp = 2.*pi*radius*um*(radius+lengthAxoncomp)*um;
      TotSurfAxon = TotSurfAxon + SurAxoncomp;
      fMassAxoncomp[fnbAxoncomp] = density*VolAxoncomp; 
      fMassAxonTot += fMassAxoncomp[fnbAxoncomp];
      Axon = Axon.unit();
   
      // Euler angles of each compartment
      G4double theta_eulerAxon = Axon.theta();
      G4double phi_eulerAxon = Axon.phi();
      G4double psi_eulerAxon = 0;

      //Rotation Matrix, Euler constructor build inverse matrix.
      G4RotationMatrix rotmAxonInv = G4RotationMatrix(phi_eulerAxon+pi/2,
						      theta_eulerAxon,
						      psi_eulerAxon);
      G4RotationMatrix rotmAxon = rotmAxonInv.inverse();
      fRotAxoncomp[fnbAxoncomp] = rotmAxon;
      ++fnbAxoncomp;
    }
    // =======================================================================
    // checking additional types
    if (typeNcomp != 1 && typeNcomp != 2 && typeNcomp != 3 && typeNcomp != 4) {
      G4cout <<  " Additional types:-->  "<< typeNcomp <<G4endl;
    }
  
    // If tracing points including spines, user can be define spine morphology
    // including stubby, mushroom, thin, long thin, filopodia and 
    // branched with heads and necks!
  
    if (typeNcomp == 5) {
      //  Sphere volume and surface area
      G4double VolSpinecomp = Piconst*pow(radius*um,3.) ;
      TotVolSpine = TotVolSpine + VolSpinecomp;
      G4double SurSpinecomp = 3.*Piconst*pow(radius*um,2.) ;
      TotSurfSpine = TotSurfSpine + SurSpinecomp;
      fMassSpinecomp[fnbSpinecomp] = density*VolSpinecomp;
      fMassSpineTot = fMassSpineTot + fMassSpinecomp[fnbSpinecomp];
      // OR    
      //  Ellipsoid volume and Approximate formula of surface area   
      // ...
      G4ThreeVector vSpine (x, y, z); 
      fPosSpinecomp[fnbSpinecomp] = vSpine; 
      fRadSpinecomp[fnbSpinecomp] = radius; 
      ++fnbSpinecomp;
    }
    ++nlines;
  }
  infile.close();
  // =======================================================================

  fnbNeuroncomp = nlines;
  G4cout <<  " Total number of compartments into Neuron : "
	 << fnbNeuroncomp << G4endl; 
  G4cout << G4endl;
  
  // to calculate SHIFT value for neuron translation
  fshiftX = (minX + maxX)/2. ;
  fshiftY = (minY + maxY)/2. ;
  fshiftZ = (minZ + maxZ)/2. ;
  
  // width, height, depth of bounding slice volume
  fwidthB  = std::fabs(minX - maxX) + maxRad;
  fheightB = std::fabs(minY - maxY) + maxRad;
  fdepthB  = std::fabs(minZ - maxZ) + maxRad;

  // diagonal length of bounding slice, that give diameter of sphere
  // for particle direction and fluence! 
  fdiagnlLength = std::sqrt(fwidthB*fwidthB + fheightB*fheightB 
			    + fdepthB*fdepthB);

  fTotVolNeuron = TotVolSoma+TotVolDend+TotVolAxon;
  fTotSurfNeuron = TotSurfSoma+TotSurfDend+TotSurfAxon;
  fTotMassNeuron = fMassSomaTot+fMassDendTot+fMassAxonTot;
 
  fTotVolSlice  = fwidthB*um*fheightB*um*fdepthB*um;
  fTotSurfSlice = 2*(fwidthB*um*fheightB*um+fheightB*um*fdepthB*um+
		     fwidthB*um*fdepthB*um);
  fTotMassSlice = 1.0 * (g/cm3) *fTotVolSlice;  

  fTotVolMedium  = Piconst*pow(fdiagnlLength*um/2.,3.) ;
  fTotSurfMedium = 3.*Piconst*pow(fdiagnlLength*um/2.,2);
  fTotMassMedium = 1.0 * (g/cm3) *fTotVolMedium;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Load prepared data file of neural network with single and multiple layers
void NeuronLoadDataFile::NeuralNetworkDATAfile(const G4String& filename)
{ 
  G4String sLine = "";
  std::ifstream infile;
  infile.open(filename.c_str());
  if (!infile.is_open()) {
    G4ExceptionDescription ed;
    ed << "Datafile " << filename << " is not opened!";
    G4Exception("NeuronLoadDataFile::NeuralNetworkDATAfile()","dna014",
		FatalException, ed, "Check file path");
    return;
  }
  G4cout << "NeuronLoadDataFile::NeuralNetworkDATAfile: opened " 
	 << filename << G4endl;  

  G4int nlines, nbSoma, nbDendrite;
  nlines = 0;
  fnbSomacomp = 0 ;  // total number of compartment into Soma 
  fnbDendritecomp = 0 ; // total number of compartment into Dendrites
  fnbAxoncomp = 0 ;  // total number of compartment into Axon
  fnbSpinecomp = 0 ; // total number of compartment into Spines 
  G4double TotVolSoma, TotVolDend, TotVolAxon;
  TotVolSoma=TotVolDend=TotVolAxon=0.;
  G4double TotSurfSoma, TotSurfDend, TotSurfAxon;
  TotSurfSoma=TotSurfDend=TotSurfAxon=0.;
  G4int typeNcomp;   // types of structure: soma, axon, apical dendrite, etc. 
  G4double x1,y1,z1,x2,y2,z2; // cartesian coordinates of each compartment 
  G4double radius; // radius of each compartment in micrometer
  G4double height; // height of each compartment in micrometer 
  G4double maxRad = -1e+09;
  G4double density = 1.0 * (g/cm3) ; // water medium
  G4double Piconst = (4.0/3.0)*pi ;
 
  for (;;) {
    getline(infile, sLine);
    if (infile.eof()) { break; }
    std::istringstream form(sLine);
    if (nlines == 0) {
      // to read total number of compartments
      form >> fnbNeuroncomp >> nbSoma >> nbDendrite ; 
      fMassSomacomp.resize(nbSoma, 0);
      fPosSomacomp.resize(nbSoma);
      fRadSomacomp.resize(nbSoma, 0); 
      fRadDendcomp.resize(nbDendrite, 0); 
      fHeightDendcomp.resize(nbDendrite, 0);
      fMassDendcomp.resize(nbDendrite, 0);
      fDistADendSoma.resize(nbDendrite, 0);
      fDistBDendSoma.resize(nbDendrite, 0);
      fPosDendcomp.resize(nbDendrite);
      fRotDendcomp.resize(nbDendrite);
    }
    // =======================================================================
    // Soma compartments represented as Sphere or Ellipsoid solid
    if (nlines > 0 && nlines <= nbSoma) {
      form >> typeNcomp >> x1 >> y1 >> z1 >> radius ;
      if (typeNcomp !=1) break;
      // max diameter of compartments 
      if (maxRad < radius) maxRad = radius;  
      //  Sphere volume and surface area
      G4double VolSomacomp = Piconst*pow(radius*um,3.) ;
      TotVolSoma = TotVolSoma + VolSomacomp;
      G4double SurSomacomp = 3.*Piconst*pow(radius*um,2.) ;
      TotSurfSoma = TotSurfSoma + SurSomacomp;
      fMassSomacomp[fnbSomacomp] = density*VolSomacomp;
      fMassSomaTot = fMassSomaTot + fMassSomacomp[fnbSomacomp];
   
      G4ThreeVector vSoma (x1 ,y1 ,z1); 
      fPosSomacomp[fnbSomacomp] = vSoma; 
      fRadSomacomp[fnbSomacomp]= radius; 
      ++fnbSomacomp; 
    }
    // =======================================================================
    // Apical and basal dendritic compartments represented as cylinderical solid    
    if (nlines > nbSoma && nlines <= fnbNeuroncomp) {
      form >> typeNcomp >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> radius >> height;
      if (typeNcomp != 3 ) break;   // || typeNcomp != 4
   
      // To calculate length, center and rotation angles of each cylinder  
      // Center-position of each cylinder
      G4double Dendxx= x1 + x2;
      G4double Dendyy= y1 + y2;
      G4double Dendzz= z1 + z2;
      G4ThreeVector translmDend = G4ThreeVector(Dendxx/2. , 
						Dendyy/2. , Dendzz/2.) ;
      fPosDendcomp [fnbDendritecomp] = translmDend;
      fRadDendcomp [fnbDendritecomp]= radius;
      G4double lengthDendcomp = height;
      // Height of compartment
      fHeightDendcomp [fnbDendritecomp]= lengthDendcomp;
      // Distance from Soma
   
      //  Cylinder volume and surface area
      G4double VolDendcomp = pi*pow(radius*um,2)*(lengthDendcomp*um);
      TotVolDend = TotVolDend + VolDendcomp;
      G4double SurDendcomp = 2.*pi*radius*um*(radius+lengthDendcomp)*um;
      TotSurfDend = TotSurfDend + SurDendcomp;
      fMassDendcomp[fnbDendritecomp] = density*VolDendcomp; 
      fMassDendTot = fMassDendTot + fMassDendcomp[fnbDendritecomp]; 
   
      G4double Dendx= x1 - x2;
      G4double Dendy= y1 - y2;
      G4double Dendz= z1 - z2; 
      Dendx=Dendx/lengthDendcomp;
      Dendy=Dendy/lengthDendcomp;
      Dendz=Dendz/lengthDendcomp;
   
      // Euler angles of each compartment
      G4ThreeVector directionDend = G4ThreeVector(Dendx,Dendy,Dendz);
      G4double theta_eulerDend = directionDend.theta();
      G4double phi_eulerDend = directionDend.phi();
      G4double psi_eulerDend = 0;

      //Rotation Matrix, Euler constructor build inverse matrix.
      G4RotationMatrix rotmDendInv = G4RotationMatrix(phi_eulerDend+pi/2,
						      theta_eulerDend,
						      psi_eulerDend);
      G4RotationMatrix rotmDend = rotmDendInv.inverse();
      
      fRotDendcomp[fnbDendritecomp] = rotmDend;
      ++fnbDendritecomp;
    }
    ++nlines;
  }
  
  // =======================================================================

  G4cout << " Total number of compartments into Neuron : " << 
    fnbNeuroncomp << G4endl;
 
  // to calculate SHIFT value for neuron translation
  fshiftX = 0.; //(minX + maxX)/2. ;
  fshiftY = 0.; //(minY + maxY)/2. ;
  fshiftZ = 0.; //(minZ + maxZ)/2. ;
 
  // width, height, depth of bounding slice volume
  fwidthB  = 640.;
  fheightB = 280.;
  fdepthB  = 25.;
  // diagonal length of bounding slice, that give diameter of sphere
  // for particle direction and fluence! 
  fdiagnlLength = std::sqrt(fwidthB*fwidthB + fheightB*fheightB 
			    + fdepthB*fdepthB);

  fTotVolNeuron = TotVolSoma+TotVolDend+TotVolAxon;
  fTotSurfNeuron = TotSurfSoma+TotSurfDend+TotSurfAxon;
  fTotMassNeuron = fMassSomaTot+fMassDendTot+fMassAxonTot;

  fTotVolSlice  = fwidthB*um*fheightB*um*fdepthB*um;
  fTotSurfSlice = 2*(fwidthB*um*fheightB*um+fheightB*um*fdepthB*um+
                  fwidthB*um*fdepthB*um);
  fTotMassSlice = 1.0 * (g/cm3) *fTotVolSlice;  

  fTotVolMedium  = Piconst*pow(fdiagnlLength*um/2.,3.) ;
  fTotSurfMedium = 3.*Piconst*pow(fdiagnlLength*um/2.,2);
  fTotMassMedium = 1.0 * (g/cm3) *fTotVolMedium; 

  infile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeuronLoadDataFile::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  // to calculate Euler angles from Rotation Matrix after Inverse!
  //
  G4RotationMatrix rotmNeuron = G4RotationMatrix(fRotNeuroncomp[copyNo]);
  G4double cosX = std::sqrt (rotmNeuron.xx()*rotmNeuron.xx() + 
			     rotmNeuron.yx()*rotmNeuron.yx()) ; 
  G4double euX, euY, euZ;
  if (cosX > 16*FLT_EPSILON) {
    euX = std::atan2 (rotmNeuron.zy(),rotmNeuron.zz());
    euY = std::atan2 (-rotmNeuron.zx(),cosX);
    euZ = std::atan2 (rotmNeuron.yx(),rotmNeuron.xx());
  } 
  else {
    euX = std::atan2 (-rotmNeuron.yz(),rotmNeuron.yy());
    euY = std::atan2 (-rotmNeuron.zx(),cosX);
    euZ = 0. ;
  }
  G4RotationMatrix* rot = new G4RotationMatrix();
  rot->rotateX(euX);
  rot->rotateY(euY);
  rot->rotateZ(euZ);  

  physVol->SetRotation(rot);  

  // shift of cylinder compartments 
  G4ThreeVector originNeuron((fPosNeuroncomp[copyNo].x()-fshiftX) * um,  
			     (fPosNeuroncomp[copyNo].y()-fshiftY) * um, 
			     (fPosNeuroncomp[copyNo].z()-fshiftZ) * um);
  physVol->SetTranslation(originNeuron);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeuronLoadDataFile::ComputeDimensions
(G4Tubs& fcylinderComp, const G4int copyNo, const G4VPhysicalVolume*) const
{ 
  fcylinderComp.SetInnerRadius(0*um);
  fcylinderComp.SetOuterRadius(fRadNeuroncomp[copyNo]*um);
  fcylinderComp.SetZHalfLength(fHeightNeuroncomp[copyNo]*um /2.);
  fcylinderComp.SetStartPhiAngle(0.*deg);
  fcylinderComp.SetDeltaPhiAngle(360.*deg); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
