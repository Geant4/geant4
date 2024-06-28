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


#ifndef NUDEXPSF_HH
#define NUDEXPSF_HH 1

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

//using namespace std;

class G4NuDEXLevelDensity;

/*
All energies in MeV
PSF are defined as in RIPL-3: PSF=Eg**(-2L-1) x Gamma width x level density
JL defines PSF x Eg**(2L+1) instead

  PSFType=0 --> SLO
  PSFType=1 --> EGLO, as defined in RIPL-3, but using always Tf in the formula
  PSFType=2 --> SMLO, as defined in RIPL-3
  PSFType=3 --> GLO  (like EGLO, but k1=k2=1)
  PSFType=4 --> MGLO (like EGLO, but k2=1)
  PSFType=5 --> KMF
  PSFType=6 --> GH
  PSFType=7 --> EGLO, but the k parameter is provided (MEGLO)
  PSFType=8 --> EGLO, but the "k1" and "k2" parameters are provided (MEGLO)
  PSFType=9 --> EGLO, but the k parameter and a constant temperature of the nucleus is provided (MEGLO)
  PSFType=10 --> EGLO, but the "k1" and "k2" parameters and a constant temperature of the nucleus are provided (MEGLO)
  PSFType=11 --> SMLO, as defined in Eur. Phys. J. A (2019) 55: 172
  PSFType=20 --> gaussian (to simulate small bumps or resonances)
  PSFType=21 --> expo -->  C*exp(-eta*Eg). It is defined with three entries: C eta dummy
  PSFType=40 --> pointwise function type 1 (only input file)
  PSFType=41 --> pointwise function type 2 (only input file)

Procedure to obtain the PSF, in order of hierarchy:
  - Get the data from inputfname
  - Get the data from PSF_param.dat file
  - Get the data from IAEA-2019 PSF values (if PSFflag==0)
  - Get the data from RIPL-3 experimental MLO values --> gdr-parameters&errors-exp-MLO.dat
  - Get the data from RIPL-3 Theorethical values --> gdr-parameters-theor.dat
  - Use RIPL-3 and RIPL-2 theoretical formulas
*/



class G4NuDEXPSF{

public:
  G4NuDEXPSF(G4int aZ,G4int aA);
  ~G4NuDEXPSF();

  G4int Init(const char* dirname,G4NuDEXLevelDensity* aLD,const char* inputfname=0,const char* defaultinputfname=0,G4int PSFflag=0);
  G4double GetE1(G4double Eg,G4double ExcitationEnergy);
  G4double GetM1(G4double Eg,G4double ExcitationEnergy);
  G4double GetE2(G4double Eg,G4double ExcitationEnergy);
  void PrintPSFParameters(std::ostream &out);
  void PrintPSFParametersInInputFileFormat(std::ostream &out);

private:

  G4bool TakePSFFromInputFile(const char* fname);
  G4bool TakePSFFromDetailedParFile(const char* fname);
  G4bool TakePSFFromIAEA01(const char* fname); // IAEA - PSF values 2019
  G4bool TakePSFFromRIPL01(const char* fname); // RIPL3-MLO values
  G4bool TakePSFFromRIPL02(const char* fname); // RIPL3-Theorethical values
  void GenerateM1AndE2FromE1(); // From RIPL-3 and RIPL-2 recommendations


  //Shapes:
  //Typical ones:
  G4double SLO(G4double Eg,G4double Er,G4double Gr,G4double sr);                          //PSFType=0
  G4double EGLO(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy); //PSFType=1
  G4double SMLO(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy); //PSFType=2
  G4double GLO(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy);  //PSFType=3
  G4double MGLO(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy); //PSFType=4
  G4double KMF(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy);  //PSFType=5
  G4double GH(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy);   //PSFType=6
  G4double MEGLO(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy,G4double k_param1,G4double k_param2,G4double Temp=-1);//PSFType=6,7,8,9,10
  G4double SMLO_v2(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy); //PSFType=11


  G4double Gauss(G4double Eg,G4double Er,G4double Gr,G4double sr); //PSFType=20
  G4double Expo(G4double Eg,G4double C,G4double eta); //PSFType=21

  //PSFType=40, PSFType=41  are pointwise defined functions
  
  //------------------------------
  G4double EGLO_GLO_MGLO(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double ExcitationEnergy,G4int Opt);
  G4double FlexibleGLOType(G4double Eg,G4double Er,G4double Gr,G4double sr,G4double Temp1,G4double k_param1,G4double Temp2,G4double k_param2);
  G4double Gamma_k(G4double Eg,G4double Er,G4double Gr,G4double Temp,G4double k_param);

private:
  G4int Z_Int,A_Int;

  G4int nR_E1,nR_M1,nR_E2;
  G4int PSFType_E1[10], PSFType_M1[10], PSFType_E2[10];
  G4double E_E1[10],G_E1[10],s_E1[10],p1_E1[10],p2_E1[10],p3_E1[10]; 
  G4double E_M1[10],G_M1[10],s_M1[10],p1_M1[10],p2_M1[10],p3_M1[10]; 
  G4double E_E2[10],G_E2[10],s_E2[10],p1_E2[10],p2_E2[10],p3_E2[10]; 

  //-----------------------------------------------
  //PSF pointwise defined PSF --> PSFType=3,4,6
  G4int np_E1,np_M1,np_E2;
  G4double *x_E1,*y_E1;
  G4double *x_M1,*y_M1;
  G4double *x_E2,*y_E2;
  G4double E1_normFac,M1_normFac,E2_normFac;
  G4double NormEmin,NormEmax;
  //-----------------------------------------------

  G4double ScaleFactor_E1,ScaleFactor_M1,ScaleFactor_E2;

  G4double EvaluateFunction(G4double xval,G4int np,G4double* x,G4double* y);
  void Renormalize();

  G4NuDEXLevelDensity* theLD;
};




#endif

