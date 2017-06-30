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
/// \file NeuronLoadDataFile.cc
/// \brief Implementation of the NeuronLoadDataFile class

#include "NeuronLoadDataFile.hh"
//#include "NeuronLoadMessenger.hh"
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
#include <cmath>
#include <sstream>
#include <string>
#include <stdlib.h>
//define if the program is running with Geant4
#define GEANT4
#ifdef GEANT4
//Specific to Geant4, globals.hh is used for G4cout
#include "globals.hh"
#endif
#include "CommandLineParser.hh"
#include "G4UImanager.hh"

using namespace std;
using namespace G4DNAPARSER;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
NeuronLoadDataFile::NeuronLoadDataFile()
{
  //CommandLineParser* parser = CommandLineParser::GetParser();
  Command* commandLine(0);
 
   // 1. Load single neuron morphology and obtain parameters.
   // Default SWC file name of neuron
   fNeuronFileNameSWC = G4String("GranuleCell-Nr2.CNG.swc");  

   // 2. Load neural network and obtain parameters.
   // Default prepared data filename of neural network with single/multi-layer.
   // Small network of 10 pyramidal neurons with single layer
   fNeuronFileNameDATA = G4String("NeuralNETWORK.dat"); 

   // Load/change SWC or DAT as "CommandLineParser" class
   if((commandLine=CommandLineParser::GetParser()->GetCommandIfActive("-swc")))
     {
       fNeuronFileNameSWC = G4String(commandLine->GetOption());
       SingleNeuronSWCfile(fNeuronFileNameSWC);  
     }
   //if (CommandLineParser::GetParser()->GetCommandIfActive("-network"))
   //  {
   //     NeuralNetworkDATAfile(fNeuronFileNameDATA); 
   //  }
   if ((commandLine =CommandLineParser::GetParser()->
                     GetCommandIfActive("-network")))
     {
       fNeuronFileNameDATA = G4String(commandLine->GetOption());
       NeuralNetworkDATAfile(fNeuronFileNameDATA); 
     }
   else 
     {     
       SingleNeuronSWCfile(fNeuronFileNameSWC);  
     }  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeuronLoadDataFile::SingleNeuronSWCfile (const G4String& filename)
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
  if (!infile)
  {
#ifdef GEANT4
    G4cout<<"\n NeuronLoadDataFile::SingleNeuronSWCfile >> datafile "
          <<filename<<" not found !!!!"<<G4endl;
    exit(0);
#endif   
  }
  else 
  {    
#ifdef G4VERBOSE 
  G4cout<<"\n NeuronLoadDataFile::SingleNeuronSWCfile >>  opening filename:  "
        << "\n" <<'\t'<<'\t'<<'\t'<<'\t'<<'\t'<<"' "<<filename 
        << " ' \n"<< G4endl;
#endif

 G4int nrows,nlines;
 nrows=0; nlines=0;
 while (getline(infile, sLine))
 {
  fnNn = new G4int[nrows];
  fpNn = new G4int[nrows];  
  fnNd = new G4int[nrows];
  fpNd = new G4int[nrows]; 
  fnNa = new G4int[nrows];
  fpNa = new G4int[nrows];   
  nrows++;
 }
 infile.close();
 //G4cout <<  " number of tracing points: "<< nrows <<G4endl;  
 infile.open(filename.c_str());
 
 fnbSomacomp = 0 ;     // total number of compartment into Soma 
 fnbDendritecomp = 0 ; // total number of compartment into Dendrites
 fnbAxoncomp = 0 ;     // total number of compartment into Axon
 fnbSpinecomp = 0 ;    // total number of compartment into Spines 
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

 fMassSomacomp  = new G4double[nrows];
 fMassSomaTot   = 0.0 ;
 fPosSomacomp   = new G4ThreeVector[nrows];
 fRadSomacomp   = new G4double[nrows]; 
 G4ThreeVector * PosDendcomp= new G4ThreeVector[nrows];
 fRadDendcomp   = new G4double[nrows]; 
 fHeightDendcomp= new G4double[nrows];
 fMassDendcomp  = new G4double[nrows];
 fMassDendTot   = 0.0 ;
 fDistADendSoma = new G4double[nrows];
 fDistBDendSoma = new G4double[nrows];
 fPosDendcomp   = new G4ThreeVector[nrows];
 fRotDendcomp   = new G4RotationMatrix[nrows];
 G4ThreeVector * PosAxoncomp= new G4ThreeVector[nrows];
 fRadAxoncomp   = new G4double[nrows]; 
 fHeightAxoncomp= new G4double[nrows];
 fMassAxoncomp  = new G4double[nrows];
 fMassAxonTot   = 0.0 ;
 fDistAxonsoma  = new G4double[nrows];
 fPosAxoncomp   = new G4ThreeVector[nrows];
 fRotAxoncomp   = new G4RotationMatrix[nrows]; 
 fMassSpinecomp = new G4double[nrows];
 fMassSpineTot  = 0.0 ;
 fPosSpinecomp   = new G4ThreeVector[nrows];
 fRadSpinecomp   = new G4double[nrows]; 
 G4ThreeVector * PosNeuroncomp= new G4ThreeVector[nrows];
 fRadNeuroncomp  = new G4double[nrows]; 
 fHeightNeuroncomp = new G4double[nrows];
 fDistNeuronsoma = new G4double[nrows];
 fPosNeuroncomp  = new G4ThreeVector[nrows];
 fRotNeuroncomp  = new G4RotationMatrix[nrows];
 fPosNeuroncomp  = new G4ThreeVector[nrows];
 fRadNeuroncomp  = new G4double[nrows]; 
 fTypeN          = new G4int[nrows];
 
 // to read datafile containing numbers, alphabets and symbols..,
 while (getline(infile, sLine))
  {
  std::istringstream form(sLine);
  G4String token;
  while (getline(form, token, ':'))
  {
   std::istringstream found(token);
   while (found >> nNcomp >> typeNcomp >> x >> y >> z >> radius >> pNcomp)
   { 
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
   if (typeNcomp == 1) 
   {   
   //  Sphere volume and surface area
   G4double VolSomacomp = Piconst*pow(radius*um,3.) ;
   TotVolSoma = TotVolSoma + VolSomacomp;
   G4double SurSomacomp = 3.*Piconst*pow(radius*um,2.) ;
   TotSurfSoma = TotSurfSoma + SurSomacomp;
   // OR    
   //  Ellipsoid volume and Approximate formula of surface area   
   //G4double VolSomacomp = Piconst*(Ra*um)*(Rb*um)*(Rc*um);   
   //G4double SurSomacomp = 3.*Piconst*pow((pow(Ra,1.6075)*pow(Rb,1.6075)+
   //pow(Ra,1.6075)*pow(Rc,1.6075)+pow(Rb,1.6075)*pow(Rc,1.6075))/3.,0.622084);
   fMassSomacomp[fnbSomacomp] = density*VolSomacomp;
   fMassSomaTot = fMassSomaTot + fMassSomacomp[fnbSomacomp];  
   G4ThreeVector vSoma (x ,y ,z); 
   fPosSomacomp [fnbSomacomp] = vSoma; 
   fRadSomacomp [fnbSomacomp]= radius; 
   // no rotate
   // OR 
   // RotationMatrix for Ellipsoid solid
   // ....
   fnbSomacomp++ ;
  } 
  // =======================================================================
  // Apical and basal dendritic compartments represented as cylinderical solid
  if (typeNcomp == 3 || typeNcomp == 4) 
  {
   G4ThreeVector vDend (x ,y ,z); 
   // Position and Radius of compartments
   PosDendcomp [fnbDendritecomp] = vDend;
   fRadDendcomp [fnbDendritecomp]= radius;
   fnNd[fnbDendritecomp]= nNcomp-(fnbSomacomp+fnbAxoncomp)-1;
   fpNd[fnbDendritecomp]= pNcomp-(fnbSomacomp+fnbAxoncomp)-1;    
   // To join two tracing points along the dendritic branches. 
   // To calculate length, center and rotation angles of each cylinder   

   // Center-position of each cylinder
   G4double Dendxx= PosDendcomp[fnNd[fnbDendritecomp]].x()+
                    PosDendcomp[fpNd[fnbDendritecomp]].x();
   G4double Dendyy= PosDendcomp[fnNd[fnbDendritecomp]].y()+
                    PosDendcomp[fpNd[fnbDendritecomp]].y();
   G4double Dendzz= PosDendcomp[fnNd[fnbDendritecomp]].z()+
                    PosDendcomp[fpNd[fnbDendritecomp]].z();
   G4ThreeVector translmDend = G4ThreeVector(Dendxx/2. , 
      Dendyy/2. , Dendzz/2.) ;
   fPosDendcomp [fnbDendritecomp] = translmDend;   
   // delta of position A and position B of cylinder 
   G4double Dendx, Dendy, Dendz;
   //primary dendritic branch should be connect with Soma
   if (fpNd[fnbDendritecomp] == -fnbSomacomp) 
   {
    Dendx= PosDendcomp[fnNd[fnbDendritecomp]].x()-
    (fPosSomacomp[0].x()+fRadSomacomp[0]);
    Dendy= PosDendcomp[fnNd[fnbDendritecomp]].y()-
    (fPosSomacomp[0].y()+fRadSomacomp[0]);
    Dendz= PosDendcomp[fnNd[fnbDendritecomp]].z()-
    (fPosSomacomp[0].z()+fRadSomacomp[0]); 
   }
   else
   {
    Dendx= PosDendcomp[fnNd[fnbDendritecomp]].x()-
           PosDendcomp[fpNd[fnbDendritecomp]].x();
    Dendy= PosDendcomp[fnNd[fnbDendritecomp]].y()-
           PosDendcomp[fpNd[fnbDendritecomp]].y();
    Dendz= PosDendcomp[fnNd[fnbDendritecomp]].z()-
           PosDendcomp[fpNd[fnbDendritecomp]].z();
   }       
   G4double lengthDendcomp = std::sqrt(Dendx*Dendx+Dendy*Dendy+Dendz*Dendz);
   // Height of compartment
   fHeightDendcomp [fnbDendritecomp]= lengthDendcomp;   
    
   // Distance from Soma
   G4double DendDisx= fPosSomacomp[0].x()-
                      fPosDendcomp [fnbDendritecomp].x();
   G4double DendDisy= fPosSomacomp[0].y()-
                      fPosDendcomp [fnbDendritecomp].y();
   G4double DendDisz= fPosSomacomp[0].z()-
                      fPosDendcomp [fnbDendritecomp].z();  
   if (typeNcomp == 3) fDistADendSoma[fnbDendritecomp] = 
      std::sqrt(DendDisx*DendDisx + DendDisy*DendDisy + DendDisz*DendDisz);
   if (typeNcomp == 4) fDistBDendSoma[fnbDendritecomp] = 
      std::sqrt(DendDisx*DendDisx + DendDisy*DendDisy + DendDisz*DendDisz);
   
   //  Cylinder volume and surface area
   G4double VolDendcomp = pi*pow(radius*um,2)*(lengthDendcomp*um);
   TotVolDend = TotVolDend + VolDendcomp;
   G4double SurDendcomp = 2.*pi*radius*um*(radius+lengthDendcomp)*um;
   TotSurfDend = TotSurfDend + SurDendcomp;
   fMassDendcomp[fnbDendritecomp] = density*VolDendcomp; 
   fMassDendTot = fMassDendTot + fMassDendcomp[fnbDendritecomp];   
   
   Dendx=Dendx/lengthDendcomp;
   Dendy=Dendy/lengthDendcomp;
   Dendz=Dendz/lengthDendcomp;
   
   // Euler angles of each compartment
   G4ThreeVector directionDend = G4ThreeVector(Dendx,Dendy,Dendz);
   G4double theta_eulerDend =  directionDend.theta();
   G4double phi_eulerDend   =  directionDend.phi();
   G4double psi_eulerDend   = 0;

   //Rotation Matrix, Euler constructor build inverse matrix.
   G4RotationMatrix rotmDendInv  = G4RotationMatrix(
      phi_eulerDend+pi/2,
      theta_eulerDend,
      psi_eulerDend);
   G4RotationMatrix rotmDend = rotmDendInv.inverse();
   /*
   // To convert from Rotation Matrix after inverse to Euler angles 
   G4double cosX = std::sqrt (rotmDend.xx()*rotmDend.xx() + 
                   rotmDend.yx()*rotmDend.yx()) ; 
   G4double euX, euY, euZ;
   if (cosX > 16*FLT_EPSILON)
    {
    euX = std::atan2 (rotmDend.zy(),rotmDend.zz());
    euY = std::atan2 (-rotmDend.zx(),cosX);
    euZ = std::atan2 (rotmDend.yx(),rotmDend.xx());
    }
   else
    {
    euX = std::atan2 (-rotmDend.yz(),rotmDend.yy());
    euY = std::atan2 (-rotmDend.zx(),cosX);
    euZ = 0. ;
    }
   G4RotationMatrix * rot = new G4RotationMatrix();
   rot->rotateX(euX);
   rot->rotateY(euY);
   rot->rotateZ(euZ);
   */
   fRotDendcomp [fnbDendritecomp]= rotmDend ;
   
   fnbDendritecomp++ ; 
  }
  
  // =======================================================================
  // Axon compartments represented as cylinderical solid
  if (typeNcomp == 2) 
  {
   G4ThreeVector vAxon (x ,y ,z); 
   // Position and Radius of compartments
   PosAxoncomp [fnbAxoncomp] = vAxon;
   fRadAxoncomp [fnbAxoncomp]= radius;
   fnNa[fnbAxoncomp]= nNcomp-(fnbSomacomp+fnbDendritecomp)-1;
   fpNa[fnbAxoncomp]= pNcomp-(fnbSomacomp+fnbDendritecomp)-1;    
   // To join two tracing points in loaded SWC data file. 
   // To calculate length, center and rotation angles of each cylinder   

   // Center-position of each cylinder
   G4double Axonxx= PosAxoncomp[fnNa[fnbAxoncomp]].x()+
                    PosAxoncomp[fpNa[fnbAxoncomp]].x();
   G4double Axonyy= PosAxoncomp[fnNa[fnbAxoncomp]].y()+
                    PosAxoncomp[fpNa[fnbAxoncomp]].y();
   G4double Axonzz= PosAxoncomp[fnNa[fnbAxoncomp]].z()+
                    PosAxoncomp[fpNa[fnbAxoncomp]].z();
   G4ThreeVector translmAxon = G4ThreeVector(Axonxx/2. , 
                               Axonyy/2. , Axonzz/2.) ;
   fPosAxoncomp [fnbAxoncomp] = translmAxon;   
   // delta of position A and position B of cylinder 
   G4double Axonx, Axony, Axonz;
   //primary axon point should be connect with Soma
   if (fpNa[fnbAxoncomp] == -(fnbSomacomp+fnbDendritecomp)) 
   {
    Axonx= PosAxoncomp[fnNa[fnbAxoncomp]].x()-
                      (fPosSomacomp[0].x()+fRadSomacomp[0]); 
    Axony= PosAxoncomp[fnNa[fnbAxoncomp]].y()-
                      (fPosSomacomp[0].y()+fRadSomacomp[0]);
    Axonz= PosAxoncomp[fnNa[fnbAxoncomp]].z()-
                      (fPosSomacomp[0].z()+fRadSomacomp[0]); 
   }
   else
   {
    Axonx= PosAxoncomp[fnNa[fnbAxoncomp]].x()-
                      PosAxoncomp[fpNa[fnbAxoncomp]].x();
    Axony= PosAxoncomp[fnNa[fnbAxoncomp]].y()-
                      PosAxoncomp[fpNa[fnbAxoncomp]].y();
    Axonz= PosAxoncomp[fnNa[fnbAxoncomp]].z()-
                      PosAxoncomp[fpNa[fnbAxoncomp]].z();
   }       
   G4double lengthAxoncomp = std::sqrt(Axonx*Axonx+Axony*Axony+Axonz*Axonz);
   // Height of compartment
   fHeightAxoncomp [fnbAxoncomp]= lengthAxoncomp;
   
   // Distance from Soma
   G4double AxonDisx= fPosSomacomp[0].x()-
                      fPosAxoncomp [fnbAxoncomp].x();
   G4double AxonDisy= fPosSomacomp[0].y()-
                      fPosAxoncomp [fnbAxoncomp].y();
   G4double AxonDisz= fPosSomacomp[0].z()-
                      fPosAxoncomp [fnbAxoncomp].z();   
   fDistAxonsoma[fnbAxoncomp] = std::sqrt(AxonDisx*AxonDisx + 
                               AxonDisy*AxonDisy + AxonDisz*AxonDisz);
         
   //  Cylinder volume and surface area
   G4double VolAxoncomp = pi*pow(radius*um,2)*(lengthAxoncomp*um);
   TotVolAxon = TotVolAxon + VolAxoncomp;
   G4double SurAxoncomp = 2.*pi*radius*um*(radius+lengthAxoncomp)*um;
   TotSurfAxon = TotSurfAxon + SurAxoncomp;
   fMassAxoncomp[fnbAxoncomp] = density*VolAxoncomp; 
   fMassAxonTot = fMassAxonTot + fMassAxoncomp[fnbAxoncomp];
   Axonx=Axonx/lengthAxoncomp;
   Axony=Axony/lengthAxoncomp;
   Axonz=Axonz/lengthAxoncomp;
   
   // Euler angles of each compartment
   G4ThreeVector directionAxon = G4ThreeVector(Axonx,Axony,Axonz);
   G4double theta_eulerAxon =  directionAxon.theta();
   G4double phi_eulerAxon   =  directionAxon.phi();
   G4double psi_eulerAxon   = 0;

   //Rotation Matrix, Euler constructor build inverse matrix.
   G4RotationMatrix rotmAxonInv  = G4RotationMatrix(
       phi_eulerAxon+pi/2,
       theta_eulerAxon,
       psi_eulerAxon);
   G4RotationMatrix rotmAxon = rotmAxonInv.inverse();
   fRotAxoncomp [fnbAxoncomp]= rotmAxon ;
   
   fnbAxoncomp++ ; 
  } 
  // =======================================================================
  // checking additional types
  if (typeNcomp != 1 && typeNcomp != 2 && typeNcomp != 3 && typeNcomp != 4)
  {
   G4cout <<  " Additional types:-->  "<< typeNcomp <<G4endl;
  }
  
  // If tracing points including spines, user can be define spine morphology
  // including stubby, mushroom, thin, long thin, filopodia and 
  // branched with heads and necks!
  
  if (typeNcomp == 5) 
  {   
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
   G4ThreeVector vSpine (x ,y ,z); 
   fPosSpinecomp [fnbSpinecomp] = vSpine; 
   fRadSpinecomp [fnbSpinecomp] = radius; 
   // no rotate
   // OR 
   // RotationMatrix for Ellipsoid solid
   // ....
   fnbSpinecomp++ ;
  }   
  
  // =======================================================================
  // Neuron- all compartments allocate in ThreeVector
   fTypeN[nlines] = typeNcomp;  
   G4ThreeVector vNeuron (x ,y ,z); 
   // Position and Radius of compartments
   PosNeuroncomp [nlines] = vNeuron;
   fRadNeuroncomp [nlines]= radius;  
   fnNn[nlines]= nNcomp-1;
   fpNn[nlines]= pNcomp-1;    
   // To join two tracing points in loaded SWC data file. 
   // To calculate length, center and rotation angles of each cylinder   

   // Center-position of each cylinder
   G4double Neuronxx= PosNeuroncomp[fnNn[nlines]].x()+
                      PosNeuroncomp[fpNn[nlines]].x();
   G4double Neuronyy= PosNeuroncomp[fnNn[nlines]].y()+
                      PosNeuroncomp[fpNn[nlines]].y();
   G4double Neuronzz= PosNeuroncomp[fnNn[nlines]].z()+
                      PosNeuroncomp[fpNn[nlines]].z();
   G4ThreeVector translmNeuron = G4ThreeVector(Neuronxx/2. , 
      Neuronyy/2. , Neuronzz/2.) ;
   fPosNeuroncomp [nlines] = translmNeuron;   
   // delta of position A and position B of cylinder 
   G4double Neuronx, Neurony, Neuronz;
   //primary point 
   if (fpNn[nlines] == -2) 
   {
    Neuronx= PosNeuroncomp[fnNn[nlines]].x()-
             fPosNeuroncomp[0].x(); 
    Neurony= PosNeuroncomp[fnNn[nlines]].y()-
             fPosNeuroncomp[0].y();
    Neuronz= PosNeuroncomp[fnNn[nlines]].z()-
             fPosNeuroncomp[0].z(); 
   }
   else
   {
    Neuronx= PosNeuroncomp[fnNn[nlines]].x()-
             PosNeuroncomp[fpNn[nlines]].x();
    Neurony= PosNeuroncomp[fnNn[nlines]].y()-
             PosNeuroncomp[fpNn[nlines]].y();
    Neuronz= PosNeuroncomp[fnNn[nlines]].z()-
             PosNeuroncomp[fpNn[nlines]].z();
   }       
   G4double lengthNeuroncomp = std::sqrt(Neuronx*Neuronx+
                               Neurony*Neurony+Neuronz*Neuronz);
   // Height of compartment
   fHeightNeuroncomp [nlines]= lengthNeuroncomp;
   // Distance from Soma
   G4double NeuronDisx= fPosNeuroncomp[0].x()-
                        fPosNeuroncomp [nlines].x();
   G4double NeuronDisy= fPosNeuroncomp[0].y()-
                        fPosNeuroncomp [nlines].y();
   G4double NeuronDisz= fPosNeuroncomp[0].z()-
                        fPosNeuroncomp [nlines].z();   
   fDistNeuronsoma[nlines] = std::sqrt(NeuronDisx*NeuronDisx + 
      NeuronDisy*NeuronDisy + NeuronDisz*NeuronDisz);
   /*      
   //  Cylinder volume and surface area
   G4double VolNeuroncomp = pi*pow(radius*um,2)*(lengthNeuroncomp*um);
   fTotVolNeuron = fTotVolNeuron + VolNeuroncomp;
   G4double SurNeuroncomp = 2.*pi*radius*um*
                            (radius+lengthNeuroncomp)*um;
   fTotSurfNeuron = TotSurfNeuron + SurNeuroncomp;
   fMassNeuroncomp[nlines] = density*VolNeuroncomp; 
   MassNeuronTot = MassNeuronTot + fMassNeuroncomp[nlines]; */
   Neuronx=Neuronx/lengthNeuroncomp;
   Neurony=Neurony/lengthNeuroncomp;
   Neuronz=Neuronz/lengthNeuroncomp;
   
   // Euler angles of each compartment
   G4ThreeVector directionNeuron = G4ThreeVector(Neuronx,Neurony,Neuronz);
   G4double theta_eulerNeuron =  directionNeuron.theta();
   G4double phi_eulerNeuron   =  directionNeuron.phi();
   G4double psi_eulerNeuron   = 0;

   //Rotation Matrix, Euler constructor build inverse matrix.
   G4RotationMatrix rotmNeuronInv  = G4RotationMatrix(
       phi_eulerNeuron+pi/2,
       theta_eulerNeuron,
       psi_eulerNeuron);
   G4RotationMatrix rotmNeuron = rotmNeuronInv.inverse();
   fRotNeuroncomp [nlines]= rotmNeuron ;   

    nlines++;
            }
        }   
    }  
    infile.close();
 // =======================================================================

 fnbNeuroncomp = nlines ;
 G4cout <<  " Total number of compartments into Neuron : "
        <<  fnbNeuroncomp<<G4endl; 
 G4cout << "\n"<<G4endl;  
  
 // to calculate SHIFT value for neuron translation
  fshiftX = (minX + maxX)/2. ;
  fshiftY = (minY + maxY)/2. ;
  fshiftZ = (minZ + maxZ)/2. ;
  
 // width, height, depth of bounding slice volume
 //maxRad = 0.0 ;
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
 
    // Soma in Violet with opacity  
    fSomaColour = new G4VisAttributes;
    fSomaColour->SetColour(G4Colour(G4Colour(0.85,0.44,0.84))); // ,1.0
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
    fNeuronColour->SetColour(G4Colour(G4Colour(0.0,0.4,0.8,0.5)));
    fNeuronColour->SetForceSolid(true);
    fNeuronColour->SetVisibility(true); 
  
   }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Load prepared data file of neural network with single and multiple layers

void NeuronLoadDataFile::NeuralNetworkDATAfile  
       (const G4String& filename)
{ 
   
  G4String sLine = "";
  std::ifstream infile;
  infile.open(filename.c_str());
  if (!infile)
  {
#ifdef GEANT4 
  G4cout<<" \n NeuronLoadDataFile::NeuralNetworkDATAfile >> datafile "
        <<filename<<" not found !!!!"<<G4endl;
        exit(0);
#endif
  }
  else 
  {    
#ifdef G4VERBOSE 
  G4cout<< " NeuronLoadDataFile::NeuralNetworkDATAfile >>  opening filename: "
        << "\n" <<'\t'<<'\t'<<'\t'<<'\t'<<'\t'<<"' "<<filename 
        << " ' \n"<< G4endl;
#endif

 G4int nlines, nbSoma, nbDendrite;
 nlines=0;
 fnbSomacomp = 0 ;  // total number of compartment into Soma 
 fnbDendritecomp = 0 ; // total number of compartment into Dendrites
 fnbAxoncomp = 0 ;  // total number of compartment into Axon
 fnbSpinecomp = 0 ; // total number of compartment into Spines 
 G4double TotVolSoma, TotVolDend, TotVolAxon;
 TotVolSoma=TotVolDend=TotVolAxon=0.;
 G4double TotSurfSoma, TotSurfDend, TotSurfAxon;
 TotSurfSoma=TotSurfDend=TotSurfAxon=0.;
 //G4int nNmorph;  // current index of neuronal morphology
 G4int typeNcomp;   // types of structure: soma, axon, apical dendrite, etc. 
 G4double x1,y1,z1,x2,y2,z2; // cartesian coordinates of each compartment 
 G4double radius; // radius of each compartment in micrometer
 G4double height; // height of each compartment in micrometer 
 //G4double minX,minY,minZ;    //minimum 
 //G4double maxX,maxY,maxZ;    //maximum 
 G4double maxRad = -1e+09;
 //minX=minY=minZ=1e+09;
 //maxX=maxY=maxZ=-1e+09;
 G4double density = 1.0 * (g/cm3) ; // water medium
 G4double Piconst = (4.0/3.0)*pi ;
 
 while (getline(infile, sLine))
 {
   std::istringstream form(sLine);
        if (nlines == 0) {
 // to read total number of compartments
 form >> fnbNeuroncomp >> nbSoma >> nbDendrite ; 
       fMassSomacomp  = new G4double[nbSoma];
 fMassSomaTot   = 0.0 ;
 fPosSomacomp   = new G4ThreeVector[nbSoma];
 fRadSomacomp   = new G4double[nbSoma]; 
 fRadDendcomp   = new G4double[nbDendrite]; 
 fHeightDendcomp = new G4double[nbDendrite];
 fMassDendcomp  = new G4double[nbDendrite];
 fMassDendTot   = 0.0 ;
 fDistADendSoma = new G4double[nbDendrite];
 fDistBDendSoma = new G4double[nbDendrite];
 fPosDendcomp   = new G4ThreeVector[nbDendrite];
 fRotDendcomp = new G4RotationMatrix[nbDendrite];
   }
   // =======================================================================
   // Soma compartments represented as Sphere or Ellipsoid solid
     if (nlines > 0 && nlines <= nbSoma) // Total number of Soma compartments
  {  
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
   // OR    
   //  Ellipsoid volume and Approximate formula of surface area   
   //G4double VolSomacomp = Piconst*(Ra*um)*(Rb*um)*(Rc*um);   
   //G4double SurSomacomp = 3.*Piconst*pow((pow(Ra,1.6075)*pow(Rb,1.6075)+
   //pow(Ra,1.6075)*pow(Rc,1.6075)+pow(Rb,1.6075)*pow(Rc,1.6075))/3.,0.622084);
   
   G4ThreeVector vSoma (x1 ,y1 ,z1); 
   fPosSomacomp [fnbSomacomp] = vSoma; 
   fRadSomacomp [fnbSomacomp]= radius; 

   // RotationMatrix for Ellipsoid solid
   // ....
   fnbSomacomp++ ; 
  }
  // =======================================================================
  // Apical and basal dendritic compartments represented as cylinderical solid    
  if (nlines > nbSoma && nlines <= fnbNeuroncomp) 
   {
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
   G4double theta_eulerDend =  directionDend.theta();
   G4double phi_eulerDend   =  directionDend.phi();
   G4double psi_eulerDend   = 0;

   //Rotation Matrix, Euler constructor build inverse matrix.
   G4RotationMatrix rotmDendInv  = G4RotationMatrix(
      phi_eulerDend+pi/2,
      theta_eulerDend,
      psi_eulerDend);
   G4RotationMatrix rotmDend = rotmDendInv.inverse();
 
   fRotDendcomp [fnbDendritecomp]= rotmDend ;
   //G4Transform3D transformDend = G4Transform3D(rotmDend,translmDend); 
   //fRotTransDendPos [fnbDendritecomp]= transformDend ;  
   fnbDendritecomp++ ; 
    
     }    
 
  nlines++;
 }
  
 // =======================================================================

 G4cout <<  " Total number of compartments into Neuron : "<< 
       fnbNeuroncomp <<G4endl; 
 G4cout << "\n"<<G4endl; 
 
 // to calculate SHIFT value for neuron translation
 fshiftX = 0.; //(minX + maxX)/2. ;
 fshiftY = 0.; //(minY + maxY)/2. ;
 fshiftZ = 0.; //(minZ + maxZ)/2. ;
 
 // width, height, depth of bounding slice volume
 //maxRad = 0.0 ;
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
 
  } 
  infile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeuronLoadDataFile::~NeuronLoadDataFile()
{
delete[] fMassSomacomp  ;
delete[] fPosSomacomp   ;
delete[] fRadSomacomp   ;
delete[] fRadDendcomp   ;
delete[] fHeightDendcomp;
delete[] fMassDendcomp  ;
delete[] fDistADendSoma ;
delete[] fDistBDendSoma ;
delete[] fPosDendcomp   ;
delete[] fRotDendcomp ;
delete[] fRadAxoncomp   ;
delete[] fHeightAxoncomp;
delete[] fMassAxoncomp  ;
delete[] fDistAxonsoma ;
delete[] fPosAxoncomp   ;
delete[] fRotAxoncomp ;
delete[] fRadNeuroncomp ;
delete[] fHeightNeuroncomp;
delete[] fMassNeuroncomp ;
delete[] fDistNeuronsoma ;
delete[] fPosNeuroncomp  ;
delete[] fRotNeuroncomp ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeuronLoadDataFile::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
/*
// sphere transformation
 G4ThreeVector
    originSoma(
      (fPosSomacomp[copyNo].x()-fshiftX) * um,  
      (fPosSomacomp[copyNo].y()-fshiftY) * um, 
      (fPosSomacomp[copyNo].z()-fshiftZ) * um 
   );
  physVol->SetRotation(0);
  physVol->SetTranslation(originSoma); 
 */ 
  
// cylinder rotation and transformation

// to calculate Euler angles from Rotation Matrix after Inverse!
//
 G4RotationMatrix rotmNeuron = G4RotationMatrix(fRotNeuroncomp[copyNo]);
 G4double cosX = std::sqrt (rotmNeuron.xx()*rotmNeuron.xx() + 
                 rotmNeuron.yx()*rotmNeuron.yx()) ; 
 G4double euX, euY, euZ;
 if (cosX > 16*FLT_EPSILON)
  {
  euX = std::atan2 (rotmNeuron.zy(),rotmNeuron.zz());
  euY = std::atan2 (-rotmNeuron.zx(),cosX);
  euZ = std::atan2 (rotmNeuron.yx(),rotmNeuron.xx());
  }
 else
  {
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
 G4ThreeVector
    originNeuron(
      (fPosNeuroncomp[copyNo].x()-fshiftX) * um,  
      (fPosNeuroncomp[copyNo].y()-fshiftY) * um, 
      (fPosNeuroncomp[copyNo].z()-fshiftZ) * um 
   );
  physVol->SetTranslation(originNeuron);  
    
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
G4VSolid* NeuronLoadDataFile::ComputeSolid(const G4int copyNo,
                                               G4VPhysicalVolume* physVol )
{
 G4VSolid* solid;
 if( typeNcomp[copyNo] == 1 ) 
 {
   solid = sphereComp ;
 }
 else if( typeNcomp[copyNo] == 3 || typeNcomp[copyNo] == 4 ) 
 {
   solid = cylinderComp ; 
 }
 return solid; 
}
*/
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
/* 
void NeuronLoadDataFile::ComputeDimensions
(G4Sphere& sphereComp, const G4int copyNo, const G4VPhysicalVolume*) const
{ 
 fsphereComp.SetInnerRadius(0);
 fsphereComp.SetOuterRadius(fRadSomacomp[copyNo] * um);
 fsphereComp.SetStartPhiAngle(0.*deg);
 fsphereComp.SetDeltaPhiAngle(360.*deg);
 fsphereComp.SetStartThetaAngle(0.*deg);
 fsphereComp.SetDeltaThetaAngle(180.*deg); 
}  
*/
/*
#if 1 
(G4Sphere& somaS, const G4int copyNo, const G4VPhysicalVolume* physVol) const
#else 
(G4Tubs& dendritesS, const G4int copyNo, const G4VPhysicalVolume* ) const
#endif 
{ 
 G4LogicalVolume* LogicalVolume = physVol->GetLogicalVolume();
 
 G4PhysicalVolume* somaPV = LogicalVolume->GetDaughter(0);
 G4LogicalVolume* somaLV = somaPV->GetLogicalVolume();
 G4Sphere* somaS = (G4Sphere*)somaLV->GetSolid();
 
 G4PhysicalVolume* dendritesPV = LogicalVolume->GetDaughter(0);
 G4LogicalVolume* dendritesLV = dendritesPV->GetLogicalVolume(); 
    G4Tubs* dendritesS = (G4Tubs*)dendritesLV->GetSolid(); 
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
