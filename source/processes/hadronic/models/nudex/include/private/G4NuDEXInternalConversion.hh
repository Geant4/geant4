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
//
// -------------------------------------------------------------------
//
//      Author:        E.Mendoza
// 
//      Creation date: May 2024
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//  NuDEX code (https://doi.org/10.1016/j.nima.2022.167894)
// 


#ifndef NUDEXINTERNALCONVERSION_HH
#define NUDEXINTERNALCONVERSION_HH 1


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

#include "G4NuDEXRandom.hh"

#define ICC_MAXNSHELLS 40
#define ICC_NMULTIP 5
#define MINZINTABLES 10 //below this value, the alpha is always 0

/*
Class to manage the internal conversion factors and the generation of converted-e-
Still not included the fluorescence-auger effects, i.e., what happens with the hole
We read the occ factors from a file, and they are stored in a matrix
The total Icc are in index=0 (data from the libraries) and index=NShells (sum of the partials)
Data are taken from: https://doi.org/10.1006/adnd.2002.0884
*/

class G4NuDEXInternalConversion{

public:
  G4NuDEXInternalConversion(G4int Z);
  ~G4NuDEXInternalConversion();
  void Init(const char* fname);
  void PrintICC(std::ostream &out);
  G4double GetICC(G4double Ene,G4int multipolarity,G4int i_shell=-1);
  G4bool SampleInternalConversion(G4double Ene,G4int multipolarity,G4double alpha=-1,G4bool CalculateProducts=true);
  void FillElectronHole(G4int i_shell); //Fluorescence/auger
  void SetRandom4Seed(unsigned int seed){theRandom4->SetSeed(seed);}


private:
  G4double Interpolate(G4double val,G4int npoints,G4double* x,G4double* y);
  void MakeTotal();


private:
  G4int theZ,NShells;
  G4double BindingEnergy[ICC_MAXNSHELLS];
  G4double *Eg[ICC_MAXNSHELLS],*Icc_E[ICC_NMULTIP][ICC_MAXNSHELLS],*Icc_M[ICC_NMULTIP][ICC_MAXNSHELLS];
  G4int np[ICC_MAXNSHELLS];
  std::string OrbitalName[ICC_MAXNSHELLS];
  G4NuDEXRandom* theRandom4;  

public:
  G4int Ne,Ng;
  G4double Eele[100],Egam[100];
};



#endif

