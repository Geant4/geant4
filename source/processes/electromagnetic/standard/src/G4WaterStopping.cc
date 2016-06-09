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
// $Id: G4WaterStopping.cc,v 1.2 2006/06/29 19:53:40 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $

//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Description: Data on stopping power
//
// Author:      V.Ivanchenko 12.05.2006
//
// Modifications:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4WaterStopping.hh"
#include "G4EmCorrections.hh"
#include "G4LPhysicsFreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4WaterStopping::G4WaterStopping(G4EmCorrections* corr)
{
  Initialise(corr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4WaterStopping::~G4WaterStopping()
{
  int n = dedx.size();
  for(int i=0; i<n; i++) {delete dedx[i];}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WaterStopping::GetElectronicDEDX(G4int iz, G4double energy)
{
  G4double res = 0.0;
  G4bool b;
  G4int idx = 0;
  for(; idx<8; idx++) {
    if(iz == Z[idx]) {
      res = (dedx[idx])->GetValue(energy, b);
      break;
    }
  }
  return res;
 }

void G4WaterStopping::Initialise(G4EmCorrections* corr)
{
  G4int i;
  //..List of ions
  G4int zz[8] = {3, 4, 5, 6, 7, 8, 9, 10 };
  G4int aa[8] = {7, 9, 11, 12, 14, 16, 19, 20 };

  for(i=0; i<8; i++) {
    Z[i] = zz[i];
    A[i] = aa[i];
  }
  //..Reduced energies
  G4double E[53] = {0.025,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,2.5,3,4,5,6,7,8,9,10,15,20,25,30,40,50,60,70,80,90,100,150,200,250
,300,400,500,600,700,800,900,1000};
  for(i=0; i<53; i++) {E[i] *= MeV;}

  //..Li
  G4double Li[53] = {2.626,2.84,3.191,3.461,3.665,3.817,3.927,4.004,4.056,4.102,3.998,3.853,3.702,3.413,3.158,2.934,2.739,2.567,2.415,2.28,1.782,1.465,1.247,1.087,0.8706,0.7299,0.631,0.5575,0.5005,0.455,0.4178,0.3004,0.2376,0.1981,0.1708,0.1354,0.1132,0.09803,0.08692,0.07842,0.07171,0.06627,0.0495,0.04085,0.03557,0.03203,0.0276,0.02498,0.02327,0.0221,0.02126,0.02064,0.02016};
  G4LPhysicsFreeVector* pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
  dedx.push_back(pv);
  for(i=0; i<53; i++) { pv->PutValues(i,E[i],Li[i]*1000.*MeV/cm); }

  //..Be
  G4double Be[53] = {3.272,3.565,4.061,4.463,4.79,5.052,5.258,5.419,5.542,5.803,5.787,5.675,5.529,5.215,4.912,4.634,4.381,4.152,3.944,3.756,3.026,2.533,2.179,1.913,1.542,1.296,1.122,0.9911,0.8898,0.8087,0.7423,0.5335,0.4219,0.3518,0.3034,0.2406,0.2013,0.1743,0.1545,0.1394,0.1275,0.1178,0.08805,0.07266,0.06328,0.05698,0.0491,0.04444,0.04141,0.03933,0.03783,0.03672,0.03588};
  pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
  dedx.push_back(pv);
  for(i=0; i<53; i++) { pv->PutValues(i,E[i],Be[i]*1000.*MeV/cm); }

  //..B
  G4double B[53] = {3.773,4.142,4.776,5.304,5.749,6.122,6.431,6.684,6.89,7.432,7.551,7.505,7.391,7.091,6.772,6.463,6.172,5.901,5.65,5.418,4.484,3.817,3.322,2.94,2.392,2.02,1.752,1.549,1.391,1.265,1.161,0.8332,0.6587,0.5492,0.4737,0.3757,0.3144,0.2723,0.2415,0.2179,0.1993,0.1842,0.1376,0.1136,0.09894,0.08909,0.07678,0.0695,0.06477,0.06151,0.05916,0.05743,0.05611};
  pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
  dedx.push_back(pv);
  for(i=0; i<53; i++) { pv->PutValues(i,E[i],B[i]*1000.*MeV/cm); }

  //..C
  G4double C[53] = {4.154,4.593,5.358,6.009,6.568,7.049,7.46,7.809,8.103,8.968,9.262,9.311,9.25,8.994,8.68,8.358,8.045,7.747,7.465,7.199,6.093,5.269,4.636,4.137,3.403,2.891,2.516,2.23,2.004,1.823,1.673,1.2,0.9483,0.7904,0.6817,0.5406,0.4525,0.392,0.3477,0.3138,0.287,0.2653,0.1983,0.1637,0.1426,0.1284,0.1107,0.1002,0.09335,0.08865,0.08528,0.08278,0.08088};
  pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
  dedx.push_back(pv);
  for(i=0; i<53; i++) { pv->PutValues(i,E[i],C[i]*1000.*MeV/cm); }

  //..N
  G4double N[53] = {4.49,4.984,5.86,6.616,7.276,7.854,8.36,8.799,9.179,10.39,10.89,11.07,11.08,10.9,10.61,10.3,9.974,9.66,9.357,9.068,7.823,6.859,6.097,5.484,4.56,3.9,3.408,3.029,2.727,2.482,2.28,1.636,1.291,1.076,0.9274,0.7354,0.6156,0.5333,0.4731,0.427,0.3906,0.3611,0.27,0.2229,0.1942,0.1749,0.1507,0.1365,0.1272,0.1208,0.1162,0.1128,0.1102};
  pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
  dedx.push_back(pv);
  for(i=0; i<53; i++) { pv->PutValues(i,E[i],N[i]*1000.*MeV/cm); }

  //..O
  G4double O[53] = {4.778,5.321,6.298,7.152,7.907,8.578,9.173,9.7,10.16,11.73,12.46,12.78,12.89,12.81,12.57,12.27,11.95,11.63,11.32,11.01,9.659,8.571,7.691,6.967,5.854,5.042,4.427,3.945,3.56,3.245,2.983,2.142,1.689,1.406,1.212,0.9602,0.8037,0.6963,0.6178,0.5577,0.5102,0.4716,0.3527,0.2913,0.2538,0.2285,0.197,0.1784,0.1663,0.1579,0.1519,0.1475,0.1441};
  pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
  dedx.push_back(pv);
  for(i=0; i<53; i++) { pv->PutValues(i,E[i],O[i]*1000.*MeV/cm); }

  //..F
  G4double F[53] = {4.992, 5.575, 6.637, 7.578, 8.418, 9.171, 9.847, 10.45, 11, 12.9, 13.88, 14.35, 14.56, 14.59, 14.4, 14.13, 13.83, 13.51, 13.19, 12.87, 11.44, 10.26, 9.279, 8.463, 7.187, 6.237, 5.506, 4.928, 4.461, 4.076, 3.753, 2.707, 2.137, 1.779, 1.533, 1.215, 1.017, 0.8809, 0.7816, 0.7056, 0.6456, 0.5969, 0.4466, 0.3688, 0.3213, 0.2894, 0.2496, 0.2259,
0.2106, 0.2, 0.1924, 0.1868, 0.1825};
  pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
  dedx.push_back(pv);
  for(i=0; i<53; i++) { pv->PutValues(i,E[i],F[i]*1000.*MeV/cm); }

  //..Ne
  G4double Ne[53] = {5.182, 5.797, 6.931, 7.948, 8.865, 9.693, 10.044, 11.12, 11.74, 13.98, 15.21, 15.85, 16.17, 16.33, 16.21, 15.98, 15.69, 15.38, 15.06, 14.74, 13.24, 11.98, 10.91, 10.01, 8.584, 7.503, 6.66, 5.986, 5.436, 4.979, 4.595, 3.332, 2.635, 2.195, 1.892, 1.499, 1.255, 1.087, 0.9646, 0.8709, 0.7969, 0.7368, 0.5514, 0.4555, 0.3969, 0.3576, 0.3083, 0.2792, 0.2602, 0.2472, 0.2378, 0.2308, 0.2255};
  pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
  dedx.push_back(pv);
  for(i=0; i<53; i++) { pv->PutValues(i,E[i],Ne[i]*1000.*MeV/cm); }

  if(corr) {
    for(i=0; i<8; i++) {corr->AddStoppingData(Z[i], A[i], "G4_WATER", *(dedx[i]));}
  }
}
