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


#ifndef NUDEXLEVELDENSITY_HH
#define NUDEXLEVELDENSITY_HH 1

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

//Level densities as they are defined in the RIPL-3 manual

//LDTYPE=1,2,3 --> Back-Shifted-Fermi-Gas model, Constant Temperature, Back-shifted: Egidy

#define DEFAULTLDTYPE 1

//using namespace std;

class G4NuDEXLevelDensity{

public:
  G4NuDEXLevelDensity(G4int aZ,G4int aA,G4int ldtype=DEFAULTLDTYPE);
  ~G4NuDEXLevelDensity(){}


  G4int ReadLDParameters(const char* dirname,const char* inputfname=0,const char* defaultinputfname=0);
  G4int CalculateLDParameters_BSFG(const char* dirname);
  G4int SearchLDParametersInInputFile(const char* inputfname);
  void GetSnD0I0Vals(G4double &aSn,G4double &aD0,G4double &aI0){aSn=Sn; aD0=D0; aI0=I0;}

  G4int GetLDType(){return LDType;}
  G4double GetNucleusTemperature(G4double ExcEnergy);
  G4double GetLevelDensity(G4double ExcEnergy_MeV,G4double spin,G4bool parity,G4bool TotalLevelDensity=false);
  G4double EstimateInverse(G4double LevDen_iMeV,G4double spin,G4bool parity); //an approximate value of ExcEnergy(rho), the inverse function of rho(ExcEnergy) - iMeV means 1/MeV
  G4double Integrate(G4double Emin,G4double Emax,G4double spin,G4bool parity);

  void PrintParameters(std::ostream &out);
  void PrintParametersInInputFileFormat(std::ostream &out);

private:

  //General info:
  G4int A_Int,Z_Int;
  G4int LDType; //=1,2,3 --> Back-Shifted-Fermi-Gas model, Constant Temperature, Back-shifted: Egidy
  G4double Sn,D0,I0; //I0 es el del nucleo A-1 (el que captura)
  G4double Ed;

  G4bool HasData;

  //Level density parameters:
  G4double A_mass,ainf_ldpar,gamma_ldpar,dW_ldpar,Delta_ldpar,T_ldpar,E0_ldpar,Ex_ldpar;


};



#endif

