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

#ifndef XLogicalLattice_h
#define XLogicalLattice_h

#define MAXRES 322  //maximum one dimensional map resolution

#include <iostream>
#include <fstream>
#include <string>
#include "G4ThreeVector.hh"


using namespace std;

class XLogicalLattice{

private:
  
  int fVresTheta; //velocity  map theta resolution (inclination)
  int fVresPhi;   //velocity  map phi resolution  (azimuth)
  int fDresTheta; //direction map theta resn
  int fDresPhi;   //direction map phi resn 

  double fMap[3][MAXRES][MAXRES];  //map for group velocity scalars
  double fN_map[3][MAXRES][MAXRES][3];
    //normalized map containing group V direction unit vectors

  double fA;       //Scaling constant for Anh.Dec. mean free path
  double fB;       //Scaling constant for Iso.Scat. mean free path
  double fDosL;    //Density of states for L-phonons
  double fDosST;   //Density of states for ST-phonons
  double fDosFT;    //Density of states for FT-phonons
  double fBeta, fGamma, fLambda, fMu; //dynamical constants for material


  ifstream fMapFile;
  int fThetaRes, fPhiRes;

public:


  
  XLogicalLattice();
  ~XLogicalLattice();

  void SetDynamicalConstants(double, double, double, double);
  void SetScatteringConstant(G4double);
  void SetAnhDecConstant(G4double);
  void SetLDOS(double);
  void SetSTDOS(double);
  void SetFTDOS(double);

    
    
  double GetBeta();
  double GetGamma();
  double GetLambda();
  double GetMu();
  G4double GetScatteringConstant();
  G4double GetAnhDecConstant();
  double GetLDOS();
  double GetSTDOS();
  double GetFTDOS();

    
  bool LoadMap(int, int, int, string);
  bool Load_NMap(int, int, int, string);
  double MapKtoV(int, G4ThreeVector);   //Get full group velocity vector
  G4ThreeVector MapKtoVDir(int, G4ThreeVector);
    //Get normalized group velocity direction so that normalisatioon
    //does not have to be done at run time

};

#endif
