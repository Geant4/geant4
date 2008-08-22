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
// $Id: G4MaterialStopping.cc,v 1.5 2008-08-22 09:23:19 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Description: Data on stopping power
//
// Author:      A.Ivantchenko 07.08.2008
//
// Modifications:
//
//----------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MaterialStopping.hh"
#include "G4EmCorrections.hh"
#include "G4LPhysicsFreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MaterialStopping::G4MaterialStopping(G4EmCorrections* corr, G4bool splineFlag) 
{
  spline = splineFlag;
  Initialise(corr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MaterialStopping::~G4MaterialStopping()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4MaterialStopping::GetDEDX(G4int ionZ, G4int idxMaterial, G4double kinEnergy)
{
  G4cout << "Hello, idxMaterial= " << idxMaterial << G4endl;
  G4double res = 0.0;
  if(ionZ < 3 || ionZ > 18 || idxMaterial < 0 || idxMaterial > 30) return res; 
  G4bool b;
  G4int idx = idxMaterial*16 + ionZ - 3;
  G4cout << "Index of DEDX " << idx << G4endl;
  G4cout << "Number  DEDX " << dedx.size() << G4endl;
  G4cout << "Pointer DEDX " << dedx[idx] << G4endl;
  G4double scaledEnergy = kinEnergy/A[ionZ - 3];
  G4cout << "Scaled Energy is " << scaledEnergy << G4endl;
  G4double emin = 0.025*MeV;
  if(scaledEnergy < emin) {
    res = (dedx[idx])->GetValue(emin, b)*sqrt(scaledEnergy/emin);
  } else {
    res = (dedx[idx])->GetValue(scaledEnergy, b);
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4MaterialStopping::GetDEDX(G4int ionZ, const G4String& NameMaterial, 
			    G4double kinEnergy)
{
  return GetDEDX(ionZ, GetMaterialIndex(NameMaterial), kinEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int
G4MaterialStopping::GetMaterialIndex(const G4String& NameMaterial)
{
  for (G4int idx=0; idx<31; idx++){
    if(MatName[idx] == NameMaterial) return idx;
  }
  return -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4MaterialStopping::GetDensity(G4int idxMaterial)
{
  if( idxMaterial < 0 || idxMaterial > 30) return 0.0;
  return Density[idxMaterial];
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String G4MaterialStopping::GetMaterialName(G4int idxMaterial)
{
  G4String s = "";
  if( idxMaterial < 0 || idxMaterial > 30) return s;
  return MatName[idxMaterial];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PhysicsVector* 
G4MaterialStopping::GetPhysicsVector(G4int ionZ, G4int idxMaterial)
{
  if(ionZ < 3 || ionZ > 18 || idxMaterial < 0 || idxMaterial > 30) return 0; 
  G4int idx = idxMaterial*16 + ionZ - 3;
  return dedx[idx];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PhysicsVector* 
G4MaterialStopping::GetPhysicsVector(G4int ionZ, const G4String& NameMaterial)
{
  return GetPhysicsVector(ionZ, GetMaterialIndex(NameMaterial));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MaterialStopping::Initialise(G4EmCorrections* corr)
{
G4int i=0, j=0;
dedx.reserve(16*31);

//..List of ions
// G4double factor = MeV*cm2/milligram;
G4double factor = 1000*MeV/cm;
G4LPhysicsFreeVector* pv;  
G4double dens[31]={1.127,.92,1.2048E-03,3.97,1.85,1.85,1.76,3.18,1.8421E-03,2.23,1.42,2.635,2.44,6.6715E-04,1.04,1.14,3.815,1.032,1.2,.94,1.4,1.19,1.06,2.20,1.8794E-03,2.32,3.667,1.0641E-03,1.8263E-03,1.0,7.5618E-04};
G4int Z_Ion[16] = {3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
G4double A_Ion[16] = {6.941,9.0122,10.811,12.011,14.007,15.999,18.998,20.180,22.990,24.305,26.982,28.086,30.974,32.065,35.453,39.948};
G4int AA_Ion[16] = {7,9,11,12,14,16,19,20,23,24,27,28,31,32,35,40};
G4String NameMaterial[31]={"G4_A-150_TISSUE","G4_ADIPOSE_TISSUE_ICRP","G4_AIR","G4_ALUMINUM_OXIDE","G4_BONE_COMPACT_ICRU","G4_BONE_CORTICAL_ICRP","G4_C-552","G4_CALCIUM_FLUORIDE","G4_CARBON_DIOXIDE","G4_Pyrex_Glass","G4_KAPTON","G4_LITHIUM_FLUORIDE","G4_LITHIUM_TETRABORATE","G4_METHANE","G4_MUSCLE_STRIATED_ICRU","G4_NYLON-6/6","G4_PHOTO_EMULSION","G4_PLASTIC_SC_VINYLTOLUENE","G4_POLYCARBONATE","G4_POLYETHYLENE","G4_MYLAR","G4_LUCITE","G4_POLYSTYRENE","G4_TEFLON","G4_PROPANE","G4_SILICON_DIOXIDE","G4_SODIUM_IODIDE","G4_TISSUE-METHANE","G4_TISSUE-PROPANE","G4_WATER","G4_WATER_VAPOR"};
//G4String Mater0[31]={"A","AA","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","XX","YY","ZZ","Z3"};

for(i=0; i<16; i++) { 
  Z[i] = Z_Ion[i];
  A[i] = A_Ion[i];
}

for(i=0;i<31;i++){
MatName[i]=NameMaterial[i];
Density[i]=dens[i]*gram/cm3;
}

G4double E[31] = {.025,.03,.04,.05,.06,.07,.08,.09,.1,.15,.2,.25,.3,.4,.5,.6,.7,.8,.9,1,1.5,2,2.5,3,4,5,6,7,8,9,10};
for(i=0;i<31;i++){E[i] *=MeV;}

// Mater0_(Z_Ion)[31]
G4double A_3[31]={2.748,2.992,3.386,3.674,3.877,4.016,4.108,4.166,4.2,4.186,4.059,3.906,3.751,3.458,3.2,2.974,2.777,2.603,2.449,2.313,1.807,1.485,1.263,1.1,.8801,.7372,.6368,.5623,.5046,.4586,.4209};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_3[i]*dens[j]*factor);}
G4double A_4[31]={3.391,3.718,4.275,4.722,5.07,5.333,5.527,5.668,5.769,5.944,5.889,5.762,5.61,5.291,4.986,4.706,4.45,4.219,4.008,3.817,3.075,2.571,2.21,1.938,1.56,1.31,1.132,1,.8972,.8151,.748};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_4[i]*dens[j]*factor);}
G4double A_5[31]={3.899,4.3,5.003,5.593,6.081,6.474,6.784,7.025,7.21,7.643,7.704,7.633,7.509,7.202,6.881,6.57,6.277,6.004,5.75,5.515,4.562,3.881,3.373,2.982,2.423,2.043,1.769,1.564,1.403,1.275,1.169};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_5[i]*dens[j]*factor);}
G4double A_6[31]={4.3,4.767,5.598,6.316,6.933,7.453,7.883,8.233,8.514,9.262,9.476,9.487,9.409,9.143,8.827,8.504,8.189,7.888,7.603,7.333,6.205,5.361,4.712,4.2,3.448,2.925,2.542,2.251,2.022,1.837,1.686};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_6[i]*dens[j]*factor);}
G4double A_7[31]={4.67,5.188,6.12,6.941,7.664,8.295,8.836,9.293,9.673,1.77,11.17,11.29,11.28,11.08,1.79,1.48,1.15,9.843,9.537,9.243,7.973,6.983,6.2,5.569,4.621,3.946,3.444,3.057,2.751,2.502,2.297};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_7[i]*dens[j]*factor);}
G4double A_8[31]={4.999,5.564,6.589,7.5,8.319,9.049,9.693,1.25,1.73,12.2,12.82,13.07,13.14,13.04,12.79,12.5,12.18,11.86,11.54,11.23,9.852,8.734,7.828,7.082,5.938,5.105,4.475,3.984,3.591,3.271,3.005};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_8[i]*dens[j]*factor);}
G4double A_9[31]={5.251,5.856,6.959,7.948,8.846,9.66,1.39,11.04,11.61,13.46,14.31,14.69,14.85,14.85,14.65,14.38,14.07,13.75,13.43,13.11,11.65,1.44,9.434,8.594,7.282,6.309,5.562,4.973,4.497,4.106,3.779};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_9[i]*dens[j]*factor);}
G4double A_10[31]={5.481,6.119,7.29,8.347,9.312,1.19,11,11.74,12.39,14.62,15.71,16.25,16.5,16.61,16.48,16.24,15.95,15.64,15.32,15,13.47,12.17,11.08,1.15,8.689,7.581,6.721,6.034,5.475,5.012,4.623};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_10[i]*dens[j]*factor);}
G4double A_11[31]={5.684,6.355,7.597,8.728,9.768,1.73,11.62,12.44,13.19,15.86,17.28,18.04,18.46,18.78,18.78,18.63,18.41,18.13,17.84,17.54,16.01,14.62,13.42,12.38,1.69,9.381,8.343,7.504,6.815,6.239,5.753};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_11[i]*dens[j]*factor);}
G4double A_12[31]={5.906,6.584,7.845,9.005,1.07,11.07,12,12.87,13.66,16.63,18.27,19.18,19.69,2.11,2.15,2.01,19.79,19.51,19.21,18.9,17.31,15.88,14.63,13.55,11.78,1.41,9.323,8.435,7.7,7.083,6.558};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_12[i]*dens[j]*factor);}
G4double A_13[31]={6.121,6.817,8.116,9.314,1.42,11.46,12.44,13.35,14.2,17.49,19.41,2.51,21.16,21.76,21.9,21.84,21.65,21.41,21.12,2.82,19.23,17.75,16.44,15.29,13.4,11.91,1.71,9.726,8.906,8.212,7.619};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_13[i]*dens[j]*factor);}
G4double A_14[31]={6.335,7.047,8.378,9.61,1.75,11.83,12.84,13.8,14.69,18.27,2.45,21.74,22.53,23.31,23.57,23.57,23.44,23.23,22.96,22.67,21.09,19.58,18.22,17.01,15,13.4,12.1,11.03,1.13,9.368,8.711};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_14[i]*dens[j]*factor);}
G4double A_15[31]={6.589,7.318,8.682,9.948,11.13,12.24,13.28,14.27,15.21,19.05,21.5,23,23.94,24.93,25.31,25.41,25.34,25.17,24.94,24.68,23.14,21.6,2.19,18.93,16.8,15.08,13.67,12.5,11.5,1.66,9.93};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_15[i]*dens[j]*factor);}
G4double A_16[31]={6.793,7.544,8.947,1.25,11.47,12.61,13.69,14.72,15.7,19.79,22.5,24.22,25.31,26.52,27.05,27.24,27.24,27.12,26.93,26.69,25.22,23.67,22.23,2.92,18.67,16.84,15.32,14.04,12.96,12.03,11.22};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_16[i]*dens[j]*factor);}
G4double A_17[31]={7.055,7.821,9.248,1.57,11.82,12.99,14.1,15.15,16.16,2.46,23.42,25.34,26.59,28.01,28.67,28.95,29.02,28.95,28.8,28.59,27.18,25.63,24.15,22.8,2.46,18.53,16.92,15.55,14.39,13.38,12.51};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_17[i]*dens[j]*factor);}
G4double A_18[31]={7.196,7.978,9.43,1.77,12.04,13.24,14.37,15.44,16.47,2.95,24.14,26.26,27.67,29.29,3.09,3.46,3.59,3.57,3.45,3.27,28.91,27.36,25.87,24.48,22.07,2.06,18.38,16.94,15.72,14.65,13.72};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],A_18[i]*dens[j]*factor);}
j++;

G4double AA_3[31]={2.852,3.097,3.49,3.778,3.983,4.125,4.221,4.283,4.321,4.313,4.182,4.022,3.859,3.553,3.284,3.05,2.845,2.665,2.507,2.366,1.846,1.516,1.289,1.123,.898,.7521,.6496,.5735,.5147,.4677,.4292};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_3[i]*dens[j]*factor);}
G4double AA_4[31]={3.53,3.861,4.421,4.867,5.214,5.479,5.677,5.823,5.929,6.119,6.064,5.93,5.769,5.434,5.115,4.822,4.557,4.317,4.1,3.902,3.139,2.623,2.254,1.977,1.591,1.336,1.155,1.019,.915,.8313,.7627};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_4[i]*dens[j]*factor);}
G4double AA_5[31]={4.062,4.473,5.185,5.776,6.262,6.655,6.968,7.213,7.405,7.862,7.928,7.854,7.721,7.396,7.057,6.731,6.426,6.142,5.879,5.635,4.655,3.957,3.438,3.04,2.469,2.082,1.804,1.594,1.431,1.3,1.192};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_5[i]*dens[j]*factor);}
G4double AA_6[31]={4.476,4.958,5.807,6.532,7.149,7.668,8.099,8.452,8.739,9.519,9.746,9.757,9.672,9.386,9.051,8.71,8.381,8.067,7.77,7.491,6.329,5.464,4.801,4.279,3.512,2.98,2.59,2.294,2.061,1.873,1.718};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_6[i]*dens[j]*factor);}
G4double AA_7[31]={4.857,5.394,6.352,7.185,7.913,8.544,9.085,9.543,9.927,11.06,11.48,11.61,11.59,11.38,11.07,1.73,1.39,1.06,9.744,9.439,8.128,7.114,6.314,5.672,4.706,4.019,3.508,3.115,2.803,2.55,2.341};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_7[i]*dens[j]*factor);}
G4double AA_8[31]={5.192,5.778,6.836,7.768,8.596,9.33,9.974,1.53,11.01,12.52,13.17,13.43,13.5,13.38,13.11,12.79,12.46,12.12,11.79,11.47,1.04,8.895,7.969,7.209,6.044,5.197,4.557,4.057,3.658,3.332,3.062};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_8[i]*dens[j]*factor);}
G4double AA_9[31]={5.447,6.073,7.214,8.232,9.145,9.967,1.7,11.35,11.92,13.81,14.69,15.09,15.25,15.24,15.02,14.73,14.4,14.06,13.73,13.39,11.88,1.63,9.606,8.75,7.413,6.423,5.664,5.064,4.581,4.183,3.85};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_9[i]*dens[j]*factor);}
G4double AA_10[31]={5.677,6.338,7.55,8.641,9.63,1.52,11.34,12.07,12.73,14.99,16.12,16.68,16.95,17.05,16.9,16.64,16.33,16,15.66,15.32,13.74,12.4,11.28,1.34,8.848,7.721,6.845,6.146,5.578,5.107,4.711};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_10[i]*dens[j]*factor);}
G4double AA_11[31]={5.885,6.579,7.866,9.036,1.1,11.08,11.99,12.81,13.56,16.27,17.73,18.53,18.96,19.28,19.27,19.1,18.85,18.55,18.24,17.92,16.33,14.9,13.67,12.6,1.87,9.542,8.487,7.634,6.934,6.349,5.854};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_11[i]*dens[j]*factor);}
G4double AA_12[31]={6.112,6.813,8.118,9.318,1.42,11.44,12.39,13.26,14.06,17.05,18.74,19.68,2.21,2.64,2.67,2.51,2.27,19.97,19.65,19.31,17.66,16.18,14.9,13.8,12,1.6,9.493,8.59,7.842,7.215,6.681};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_12[i]*dens[j]*factor);}
G4double AA_13[31]={6.332,7.05,8.393,9.633,1.78,11.85,12.84,13.77,14.63,17.94,19.9,21.05,21.72,22.33,22.47,22.38,22.18,21.91,21.61,21.29,19.63,18.1,16.75,15.58,13.64,12.12,1.9,9.904,9.069,8.364,7.761};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_13[i]*dens[j]*factor);}
G4double AA_14[31]={6.551,7.286,8.66,9.935,11.12,12.23,13.26,14.23,15.14,18.74,2.97,22.31,23.13,23.93,24.18,24.17,24.02,23.78,23.5,23.19,21.53,19.96,18.57,17.33,15.28,13.65,12.33,11.23,1.32,9.542,8.873};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_14[i]*dens[j]*factor);}
G4double AA_15[31]={6.812,7.564,8.97,1.27,11.5,12.64,13.72,14.73,15.68,19.54,22.04,23.59,24.57,25.59,25.97,26.05,25.96,25.77,25.52,25.23,23.62,22.02,2.58,19.29,17.11,15.35,13.92,12.72,11.71,1.85,1.11};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_15[i]*dens[j]*factor);}
G4double AA_16[31]={7.022,7.795,9.24,1.58,11.84,13.02,14.14,15.19,16.19,2.31,23.06,24.83,25.98,27.22,27.75,27.92,27.9,27.77,27.55,27.3,25.73,24.13,22.64,21.3,19.01,17.14,15.59,14.29,13.19,12.25,11.43};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_16[i]*dens[j]*factor);}
G4double AA_17[31]={7.291,8.078,9.546,1.91,12.2,13.41,14.55,15.64,16.67,21,24,25.98,27.28,28.74,29.41,29.68,29.73,29.64,29.47,29.24,27.74,26.13,24.61,23.22,2.83,18.86,17.22,15.83,14.65,13.62,12.73};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_17[i]*dens[j]*factor);}
G4double AA_18[31]={7.435,8.239,9.732,11.11,12.42,13.66,14.83,15.94,16.99,21.52,24.74,26.92,28.38,3.07,3.87,31.24,31.35,31.31,31.17,3.97,29.52,27.91,26.37,24.95,22.48,2.43,18.71,17.25,16,14.92,13.97};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],AA_18[i]*dens[j]*factor);}
j++;

G4double B_3[31]={2.062,2.238,2.532,2.764,2.943,3.077,3.176,3.247,3.296,3.354,3.282,3.176,3.061,2.839,2.639,2.461,2.305,2.167,2.044,1.933,1.523,1.257,1.072,.9372,.752,.6313,.5463,.4831,.4341,.3949,.3627};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_3[i]*dens[j]*factor);}
G4double B_4[31]={2.56,2.797,3.205,3.546,3.829,4.06,4.245,4.39,4.502,4.751,4.76,4.686,4.581,4.346,4.113,3.896,3.696,3.513,3.346,3.193,2.592,2.179,1.879,1.653,1.335,1.123,.9729,.8601,.7724,.7024,.645};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_4[i]*dens[j]*factor);}
G4double B_5[31]={2.959,3.253,3.766,4.204,4.582,4.905,5.178,5.404,5.588,6.086,6.216,6.203,6.13,5.915,5.676,5.439,5.212,4.998,4.799,4.612,3.848,3.292,2.872,2.547,2.076,1.755,1.522,1.347,1.21,1.1,1.009};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_5[i]*dens[j]*factor);}
G4double B_6[31]={3.28,3.628,4.24,4.769,5.234,5.642,5.998,6.305,6.566,7.35,7.635,7.707,7.683,7.514,7.286,7.045,6.805,6.572,6.35,6.138,5.239,4.553,4.019,3.593,2.962,2.518,2.192,1.943,1.746,1.588,1.458};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_6[i]*dens[j]*factor);}
G4double B_7[31]={3.554,3.946,4.645,5.254,5.796,6.28,6.712,7.094,7.428,8.517,8.989,9.172,9.216,9.119,8.923,8.692,8.451,8.21,7.974,7.746,6.741,5.941,5.298,4.775,3.979,3.407,2.978,2.646,2.382,2.168,1.99};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_7[i]*dens[j]*factor);}
G4double B_8[31]={3.796,4.228,5.005,5.69,6.302,6.856,7.357,7.808,8.212,9.609,1.28,1.59,1.72,1.72,1.57,1.36,1.13,9.893,9.653,9.415,8.332,7.434,6.695,6.079,5.121,4.415,3.877,3.456,3.117,2.84,2.61};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_8[i]*dens[j]*factor);}
G4double B_9[31]={4.002,4.467,5.312,6.064,6.739,7.354,7.915,8.427,8.892,1.58,11.47,11.91,12.13,12.23,12.13,11.95,11.73,11.5,11.25,11.01,9.875,8.901,8.082,7.389,6.291,5.466,4.828,4.321,3.91,3.572,3.288};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_9[i]*dens[j]*factor);}
G4double B_10[31]={4.19,4.683,5.587,6.399,7.133,7.803,8.418,8.986,9.507,11.47,12.58,13.18,13.49,13.71,13.67,13.52,13.32,13.1,12.86,12.61,11.43,1.39,9.506,8.744,7.517,6.578,5.842,5.252,4.768,4.367,4.029};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_10[i]*dens[j]*factor);}
G4double B_11[31]={4.337,4.854,5.813,6.684,7.477,8.206,8.878,9.502,1.08,12.34,13.71,14.5,14.94,15.33,15.4,15.33,15.18,14.98,14.77,14.54,13.38,12.3,11.34,1.51,9.131,8.049,7.183,6.477,5.894,5.404,4.988};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_11[i]*dens[j]*factor);}
G4double B_12[31]={4.499,5.029,6.02,6.929,7.763,8.531,9.243,9.905,1.52,13,14.57,15.51,16.07,16.58,16.7,16.66,16.52,16.33,16.11,15.87,14.66,13.53,12.52,11.64,1.18,9.031,8.103,7.343,6.711,6.178,5.723};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_12[i]*dens[j]*factor);}
G4double B_13[31]={4.647,5.193,6.215,7.16,8.035,8.846,9.599,1.3,1.96,13.65,15.45,16.57,17.25,17.93,18.16,18.17,18.08,17.91,17.71,17.49,16.28,15.12,14.07,13.14,11.57,1.32,9.31,8.47,7.765,7.166,6.653};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_13[i]*dens[j]*factor);}
G4double B_14[31]={4.794,5.352,6.401,7.378,8.289,9.137,9.928,1.66,11.36,14.25,16.25,17.55,18.37,19.22,19.55,19.64,19.58,19.45,19.27,19.06,17.87,16.68,15.6,14.62,12.96,11.62,1.52,9.613,8.84,8.18,7.611};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_14[i]*dens[j]*factor);}
G4double B_15[31]={4.946,5.515,6.588,7.595,8.539,9.422,1.24,11.02,11.75,14.83,17.04,18.51,19.48,2.53,2.99,21.16,21.17,21.08,2.93,2.75,19.6,18.41,17.29,16.28,14.53,13.09,11.9,1.9,1.05,9.324,8.69};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_15[i]*dens[j]*factor);}
G4double B_16[31]={5.082,5.665,6.766,7.802,8.779,9.697,1.55,11.36,12.13,15.38,17.79,19.45,2.56,21.82,22.41,22.67,22.74,22.7,22.59,22.43,21.36,2.17,19.03,17.98,16.15,14.62,13.34,12.26,11.33,1.53,9.839};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_16[i]*dens[j]*factor);}
G4double B_17[31]={5.242,5.836,6.957,8.015,9.017,9.963,1.85,11.69,12.48,15.89,18.48,2.31,21.57,23.04,23.77,24.11,24.25,24.25,24.18,24.05,23.04,21.86,2.7,19.62,17.71,16.11,14.75,13.6,12.6,11.73,1.98};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_17[i]*dens[j]*factor);}
G4double B_18[31]={5.351,5.954,7.092,8.167,9.189,1.15,11.07,11.93,12.75,16.29,19.03,21.03,22.44,24.11,24.96,25.39,25.59,25.64,25.59,25.49,24.53,23.35,22.18,21.08,19.11,17.45,16.03,14.81,13.76,12.85,12.04};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],B_18[i]*dens[j]*factor);}
j++;

G4double C_3[31]={1.355,1.466,1.653,1.806,1.931,2.032,2.115,2.181,2.234,2.369,2.391,2.366,2.32,2.204,2.082,1.966,1.859,1.76,1.671,1.59,1.276,1.066,.9168,.8057,.6519,.5505,.4785,.4246,.3826,.3489,.3212};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_3[i]*dens[j]*factor);}
G4double C_4[31]={1.73,1.882,2.144,2.364,2.553,2.714,2.851,2.967,3.065,3.355,3.462,3.485,3.466,3.368,3.24,3.105,2.973,2.847,2.728,2.617,2.166,1.843,1.603,1.419,1.156,.979,.8516,.7556,.6807,.6205,.5712};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_4[i]*dens[j]*factor);}
G4double C_5[31]={2.04,2.23,2.562,2.845,3.093,3.311,3.503,3.671,3.816,4.292,4.511,4.603,4.627,4.576,4.464,4.328,4.185,4.043,3.905,3.773,3.207,2.777,2.444,2.182,1.796,1.528,1.332,1.183,1.066,.9714,.8937};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_5[i]*dens[j]*factor);}
G4double C_6[31]={2.295,2.523,2.922,3.267,3.572,3.844,4.087,4.306,4.5,5.177,5.527,5.704,5.786,5.803,5.722,5.599,5.458,5.31,5.161,5.015,4.36,3.835,3.415,3.074,2.559,2.191,1.917,1.706,1.538,1.403,1.29};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_6[i]*dens[j]*factor);}
G4double C_7[31]={2.515,2.773,3.232,3.633,3.99,4.312,4.604,4.869,5.109,5.992,6.491,6.769,6.923,7.029,6.998,6.901,6.771,6.626,6.475,6.322,5.603,4.998,4.497,4.081,3.436,2.963,2.603,2.323,2.098,1.915,1.762};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_7[i]*dens[j]*factor);}
G4double C_8[31]={2.707,2.993,3.506,3.959,4.366,4.735,5.073,5.383,5.667,6.753,7.409,7.8,8.035,8.247,8.278,8.218,8.11,7.976,7.829,7.675,6.915,6.244,5.673,5.188,4.416,3.836,3.388,3.033,2.745,2.509,2.311};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_8[i]*dens[j]*factor);}
G4double C_9[31]={2.876,3.185,3.746,4.245,4.696,5.108,5.487,5.838,6.162,7.438,8.249,8.755,9.075,9.399,9.497,9.478,9.395,9.276,9.136,8.985,8.201,7.48,6.852,6.309,5.429,4.753,4.222,3.795,3.447,3.157,2.914};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_9[i]*dens[j]*factor);}
G4double C_10[31]={3.029,3.359,3.961,4.502,4.994,5.446,5.863,6.251,6.613,8.069,9.037,9.665,1.07,1.52,1.69,1.72,1.67,1.57,1.44,1.3,9.503,8.742,8.065,7.471,6.491,5.724,5.113,4.617,4.207,3.864,3.574};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_10[i]*dens[j]*factor);}
G4double C_11[31]={3.151,3.498,4.137,4.716,5.247,5.736,6.191,6.615,7.012,8.647,9.778,1.54,11.05,11.63,11.89,11.98,11.97,11.9,11.78,11.65,1.86,1.07,9.354,8.713,7.637,6.781,6.089,5.52,5.046,4.647,4.306};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_11[i]*dens[j]*factor);}
G4double C_12[31]={3.282,3.642,4.308,4.919,5.482,6.005,6.491,6.947,7.375,9.168,1.45,11.34,11.95,12.68,13.04,13.19,13.22,13.17,13.08,12.96,12.18,11.37,1.62,9.939,8.782,7.848,7.083,6.449,5.916,5.463,5.073};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_12[i]*dens[j]*factor);}
G4double C_13[31]={3.4,3.772,4.463,5.1,5.691,6.244,6.762,7.248,7.705,9.646,11.08,12.1,12.81,13.7,14.16,14.38,14.46,14.45,14.38,14.28,13.53,12.71,11.93,11.21,9.982,8.972,8.136,7.436,6.844,6.336,5.898};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_13[i]*dens[j]*factor);}
G4double C_14[31]={3.513,3.896,4.608,5.268,5.884,6.463,7.008,7.521,8.006,1.08,11.65,12.8,13.62,14.66,15.23,15.53,15.66,15.69,15.65,15.57,14.85,14.02,13.22,12.48,11.18,1.1,9.201,8.441,7.792,7.234,6.748};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_14[i]*dens[j]*factor);}
G4double C_15[31]={3.621,4.013,4.743,5.423,6.063,6.666,7.237,7.776,8.287,1.49,12.2,13.48,14.41,15.63,16.32,16.71,16.91,16.99,16.98,16.93,16.29,15.47,14.65,13.88,12.52,11.37,1.4,9.571,8.86,8.243,7.705};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_15[i]*dens[j]*factor);}
G4double C_16[31]={3.72,4.122,4.874,5.575,6.238,6.865,7.46,8.024,8.561,1.89,12.73,14.14,15.19,16.58,17.39,17.87,18.15,18.28,18.32,18.29,17.73,16.93,16.11,15.32,13.9,12.69,11.65,1.76,9.984,9.311,8.719};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_16[i]*dens[j]*factor);}
G4double C_17[31]={3.832,4.242,5.008,5.726,6.406,7.054,7.669,8.255,8.814,11.26,13.22,14.75,15.91,17.47,18.41,18.99,19.33,19.52,19.6,19.61,19.13,18.35,17.52,16.72,15.25,13.98,12.88,11.93,11.1,1.38,9.736};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_17[i]*dens[j]*factor);}
G4double C_18[31]={3.921,4.338,5.118,5.849,6.544,7.208,7.842,8.448,9.026,11.58,13.65,15.29,16.55,18.29,19.35,2.02,2.43,2.67,2.8,2.84,2.43,19.67,18.84,18.02,16.51,15.18,14.03,13.04,12.16,11.39,1.71};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],C_18[i]*dens[j]*factor);}
j++;

G4double D_3[31]={2.204,2.403,2.723,2.962,3.137,3.265,3.355,3.418,3.459,3.491,3.406,3.291,3.171,2.939,2.731,2.548,2.387,2.244,2.117,2.004,1.581,1.307,1.117,.9771,.7854,.6601,.5718,.5059,.4547,.4138,.3803};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_3[i]*dens[j]*factor);}
G4double D_4[31]={2.712,2.981,3.438,3.805,4.096,4.324,4.5,4.636,4.738,4.953,4.942,4.856,4.742,4.494,4.252,4.026,3.82,3.631,3.459,3.301,2.684,2.259,1.951,1.719,1.391,1.172,1.016,.8994,.8084,.7355,.6757};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_4[i]*dens[j]*factor);}
G4double D_5[31]={3.108,3.438,4.018,4.502,4.904,5.236,5.506,5.725,5.9,6.355,6.459,6.43,6.345,6.115,5.863,5.616,5.381,5.16,4.954,4.762,3.975,3.404,2.974,2.64,2.158,1.827,1.588,1.406,1.264,1.15,1.057};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_5[i]*dens[j]*factor);}
G4double D_6[31]={3.419,3.802,4.489,5.081,5.589,6.02,6.385,6.691,6.946,7.683,7.935,7.988,7.949,7.761,7.52,7.267,7.017,6.776,6.546,6.328,5.402,4.698,4.151,3.715,3.069,2.615,2.28,2.024,1.822,1.658,1.524};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_6[i]*dens[j]*factor);}
G4double D_7[31]={3.699,4.122,4.892,5.573,6.171,6.693,7.146,7.536,7.87,8.913,9.345,9.503,9.529,9.41,9.198,8.955,8.703,8.452,8.208,7.972,6.937,6.116,5.46,4.925,4.113,3.528,3.09,2.75,2.479,2.259,2.077};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_7[i]*dens[j]*factor);}
G4double D_8[31]={3.948,4.407,5.25,6.009,6.688,7.293,7.829,8.3,8.711,1.07,1.7,1.99,11.09,11.06,1.9,1.67,1.43,1.18,9.928,9.682,8.563,7.642,6.886,6.257,5.281,4.562,4.014,3.583,3.237,2.954,2.718};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_8[i]*dens[j]*factor);}
G4double D_9[31]={4.148,4.637,5.539,6.363,7.111,7.788,8.396,8.94,9.423,11.09,11.94,12.35,12.54,12.61,12.5,12.3,12.07,11.82,11.57,11.32,1.15,9.147,8.308,7.602,6.483,5.644,4.994,4.477,4.058,3.711,3.421};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_9[i]*dens[j]*factor);}
G4double D_10[31]={4.332,4.847,5.801,6.677,7.485,8.224,8.897,9.508,1.06,12.03,13.09,13.65,13.94,14.13,14.07,13.91,13.7,13.46,13.22,12.96,11.75,1.68,9.773,8.995,7.745,6.79,6.041,5.439,4.946,4.535,4.189};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_10[i]*dens[j]*factor);}
G4double D_11[31]={4.494,5.033,6.04,6.97,7.838,8.642,9.383,1.06,1.69,13,14.34,15.1,15.53,15.9,15.95,15.87,15.71,15.5,15.28,15.04,13.83,12.71,11.72,1.86,9.432,8.319,7.429,6.705,6.106,5.604,5.177};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_11[i]*dens[j]*factor);}
G4double D_12[31]={4.672,5.218,6.244,7.198,8.094,8.932,9.712,1.43,11.1,13.66,15.2,16.1,16.62,17.1,17.21,17.15,17,16.8,16.58,16.34,15.09,13.92,12.89,11.99,1.5,9.319,8.373,7.597,6.952,6.407,5.942};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_12[i]*dens[j]*factor);}
G4double D_13[31]={4.842,5.403,6.458,7.44,8.367,9.243,1.06,1.83,11.55,14.36,16.13,17.2,17.86,18.51,18.72,18.72,18.62,18.44,18.24,18.01,16.76,15.57,14.49,13.53,11.93,1.66,9.617,8.759,8.039,7.428,6.903};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_13[i]*dens[j]*factor);}
G4double D_14[31]={5.012,5.585,6.664,7.673,8.626,9.531,1.39,11.19,11.95,14.99,16.98,18.23,19.02,19.84,2.16,2.23,2.17,2.03,19.84,19.63,18.4,17.18,16.06,15.06,13.37,12,1.88,9.941,9.152,8.477,7.895};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_14[i]*dens[j]*factor);}
G4double D_15[31]={5.2,5.786,6.889,7.923,8.9,9.832,1.72,11.56,12.35,15.61,17.82,19.26,2.19,21.2,21.64,21.8,21.79,21.7,21.54,21.35,2.17,18.94,17.8,16.75,14.96,13.5,12.28,11.26,1.39,9.648,9.002};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_15[i]*dens[j]*factor);}
G4double D_16[31]={5.36,5.962,7.094,8.155,9.159,1.12,11.03,11.91,12.73,16.2,18.62,2.24,21.32,22.54,23.11,23.36,23.42,23.37,23.25,23.09,21.97,2.74,19.57,18.49,16.61,15.06,13.75,12.65,11.7,1.89,1.18};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_16[i]*dens[j]*factor);}
G4double D_17[31]={5.554,6.169,7.319,8.4,9.423,1.4,11.34,12.24,13.09,16.74,19.36,21.16,22.38,23.8,24.5,24.83,24.96,24.96,24.88,24.74,23.69,22.47,21.28,2.17,18.22,16.58,15.19,14.02,13,12.12,11.35};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_17[i]*dens[j]*factor);}
G4double D_18[31]={5.672,6.299,7.469,8.566,9.605,1.6,11.55,12.47,13.35,17.15,19.95,21.92,23.28,24.91,25.73,26.16,26.34,26.39,26.34,26.22,25.23,24.02,22.82,21.68,19.67,17.97,16.52,15.28,14.21,13.28,12.46};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],D_18[i]*dens[j]*factor);}
j++;

G4double E_3[31]={2.041,2.227,2.528,2.754,2.923,3.048,3.138,3.201,3.244,3.289,3.217,3.113,3.003,2.789,2.596,2.425,2.275,2.141,2.021,1.915,1.515,1.256,1.074,.9407,.7572,.6371,.5522,.4889,.4397,.4003,.3679};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_3[i]*dens[j]*factor);}
G4double E_4[31]={2.507,2.76,3.189,3.535,3.811,4.031,4.203,4.337,4.44,4.665,4.667,4.593,4.491,4.264,4.041,3.831,3.639,3.463,3.301,3.153,2.572,2.169,1.876,1.654,1.341,1.132,.9817,.8693,.7817,.7114,.6538};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_4[i]*dens[j]*factor);}
G4double E_5[31]={2.869,3.179,3.723,4.178,4.559,4.875,5.135,5.349,5.521,5.982,6.098,6.081,6.009,5.801,5.571,5.343,5.125,4.92,4.727,4.547,3.808,3.268,2.859,2.541,2.08,1.764,1.534,1.359,1.223,1.113,1.023};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_5[i]*dens[j]*factor);}
G4double E_6[31]={3.155,3.512,4.157,4.714,5.192,5.601,5.95,6.245,6.492,7.226,7.489,7.554,7.528,7.364,7.145,6.914,6.683,6.46,6.246,6.042,5.174,4.51,3.991,3.576,2.959,2.524,2.203,1.957,1.762,1.605,1.475};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_6[i]*dens[j]*factor);}
G4double E_7[31]={3.408,3.802,4.524,5.165,5.729,6.222,6.653,7.026,7.348,8.376,8.816,8.984,9.023,8.928,8.74,8.519,8.288,8.057,7.831,7.612,6.644,5.871,5.249,4.741,3.966,3.407,2.986,2.659,2.399,2.186,2.011};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_7[i]*dens[j]*factor);}
G4double E_8[31]={3.634,4.062,4.85,5.564,6.204,6.775,7.283,7.731,8.126,9.455,1.09,1.38,1.5,1.5,1.35,1.15,9.932,9.702,9.471,9.243,8.2,7.334,6.619,6.023,5.092,4.405,3.879,3.465,3.133,2.859,2.632};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_8[i]*dens[j]*factor);}
G4double E_9[31]={3.82,4.275,5.117,5.891,6.596,7.235,7.81,8.326,8.787,1.41,11.25,11.67,11.88,11.97,11.88,11.71,11.5,11.27,11.04,1.81,9.719,8.781,7.99,7.319,6.254,5.451,4.828,4.331,3.928,3.594,3.314};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_9[i]*dens[j]*factor);}
G4double E_10[31]={3.991,4.47,5.359,6.182,6.943,7.64,8.277,8.855,9.378,11.29,12.34,12.91,13.21,13.41,13.38,13.24,13.06,12.84,12.62,12.38,11.26,1.26,9.4,8.663,7.473,6.56,5.842,5.263,4.789,4.393,4.059};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_10[i]*dens[j]*factor);}
G4double E_11[31]={4.141,4.641,5.578,6.448,7.264,8.022,8.721,9.364,9.953,12.18,13.49,14.24,14.68,15.06,15.14,15.07,14.93,14.75,14.55,14.33,13.21,12.17,11.24,1.42,9.075,8.016,7.167,6.475,5.901,5.419,5.009};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_11[i]*dens[j]*factor);}
G4double E_12[31]={4.304,4.812,5.769,6.662,7.505,8.297,9.034,9.718,1.35,12.8,14.31,15.21,15.74,16.24,16.37,16.33,16.2,16.03,15.83,15.61,14.45,13.37,12.4,11.54,1.13,9.002,8.096,7.352,6.731,6.207,5.758};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_12[i]*dens[j]*factor);}
G4double E_13[31]={4.461,4.982,5.965,6.884,7.756,8.583,9.36,1.09,1.76,13.45,15.18,16.25,16.91,17.57,17.8,17.83,17.74,17.59,17.41,17.2,16.06,14.94,13.93,13.03,11.51,1.29,9.299,8.476,7.784,7.196,6.69};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_13[i]*dens[j]*factor);}
G4double E_14[31]={4.617,5.149,6.155,7.098,7.993,8.848,9.658,1.42,11.14,14.04,15.98,17.22,18.01,18.84,19.17,19.26,19.23,19.11,18.94,18.75,17.63,16.49,15.45,14.5,12.89,11.59,1.52,9.621,8.862,8.214,7.653};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_14[i]*dens[j]*factor);}
G4double E_15[31]={4.786,5.33,6.357,7.322,8.239,9.119,9.959,1.75,11.51,14.61,16.76,18.17,19.1,2.13,2.58,2.76,2.77,2.7,2.57,2.4,19.33,18.19,17.11,16.13,14.43,13.04,11.88,1.9,1.07,9.349,8.727};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_15[i]*dens[j]*factor);}
G4double E_16[31]={4.931,5.49,6.543,7.534,8.475,9.379,1.25,11.07,11.86,15.16,17.5,19.1,2.17,21.39,21.98,22.24,22.32,22.29,22.2,22.05,21.04,19.91,18.82,17.8,16.03,14.54,13.3,12.24,11.34,1.55,9.865};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_16[i]*dens[j]*factor);}
G4double E_17[31]={5.106,5.676,6.746,7.755,8.714,9.635,1.52,11.37,12.19,15.65,18.18,19.95,21.16,22.59,23.3,23.65,23.79,23.81,23.75,23.63,22.69,21.57,2.46,19.42,17.57,16.02,14.7,13.57,12.6,11.75,11.01};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_17[i]*dens[j]*factor);}
G4double E_18[31]={5.216,5.797,6.886,7.91,8.884,9.818,1.72,11.59,12.43,16.04,18.74,2.67,22.02,23.64,24.47,24.91,25.12,25.18,25.15,25.06,24.18,23.06,21.94,2.88,18.98,17.36,15.98,14.8,13.77,12.87,12.09};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],E_18[i]*dens[j]*factor);}
j++;

G4double F_3[31]={1.926,2.103,2.397,2.623,2.794,2.92,3.013,3.079,3.126,3.192,3.136,3.044,2.941,2.737,2.551,2.384,2.236,2.105,1.987,1.882,1.486,1.229,1.05,.9178,.7373,.6196,.5365,.4747,.4267,.3883,.3569};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_3[i]*dens[j]*factor);}
G4double F_4[31]={2.371,2.606,3.017,3.357,3.634,3.855,4.029,4.165,4.271,4.519,4.545,4.488,4.397,4.186,3.972,3.77,3.582,3.409,3.251,3.105,2.528,2.129,1.838,1.618,1.309,1.102,.955,.8448,.7592,.6906,.6344};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_4[i]*dens[j]*factor);}
G4double F_5[31]={2.721,3.009,3.521,3.961,4.337,4.653,4.915,5.129,5.304,5.787,5.931,5.937,5.881,5.695,5.478,5.26,5.048,4.848,4.659,4.482,3.751,3.214,2.808,2.492,2.034,1.721,1.494,1.323,1.189,1.081,.9929};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_5[i]*dens[j]*factor);}
G4double F_6[31]={2.999,3.334,3.936,4.466,4.933,5.339,5.687,5.983,6.233,6.986,7.28,7.372,7.366,7.23,7.029,6.809,6.587,6.37,6.161,5.961,5.103,4.443,3.926,3.513,2.899,2.468,2.15,1.907,1.715,1.561,1.433};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_6[i]*dens[j]*factor);}
G4double F_7[31]={3.243,3.614,4.29,4.894,5.437,5.922,6.35,6.724,7.048,8.093,8.565,8.766,8.83,8.769,8.602,8.395,8.174,7.951,7.731,7.516,6.56,5.791,5.171,4.665,3.893,3.336,2.919,2.595,2.338,2.129,1.956};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_7[i]*dens[j]*factor);}
G4double F_8[31]={3.463,3.867,4.607,5.277,5.887,6.443,6.944,7.391,7.787,9.131,9.799,1.13,1.27,1.31,1.19,1.01,9.799,9.579,9.355,9.133,8.105,7.245,6.531,5.936,5.007,4.322,3.798,3.388,3.058,2.788,2.563};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_8[i]*dens[j]*factor);}
G4double F_9[31]={3.645,4.075,4.87,5.595,6.263,6.877,7.441,7.953,8.413,1.05,1.91,11.37,11.61,11.75,11.68,11.53,11.34,11.12,1.9,1.68,9.598,8.667,7.878,7.21,6.147,5.347,4.727,4.234,3.834,3.504,3.228};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_9[i]*dens[j]*factor);}
G4double F_10[31]={3.812,4.265,5.107,5.881,6.597,7.263,7.881,8.45,8.97,1.88,11.96,12.56,12.9,13.15,13.15,13.04,12.86,12.66,12.44,12.22,11.11,1.11,9.261,8.527,7.341,6.432,5.718,5.144,4.674,4.283,3.953};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_10[i]*dens[j]*factor);}
G4double F_11[31]={3.948,4.421,5.31,6.132,6.899,7.618,8.292,8.921,9.504,11.73,13.06,13.85,14.32,14.76,14.87,14.83,14.71,14.54,14.35,14.14,13.05,12.02,11.1,1.29,8.949,7.893,7.048,6.358,5.787,5.308,4.901};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_11[i]*dens[j]*factor);}
G4double F_12[31]={4.1,4.581,5.489,6.337,7.131,7.879,8.584,9.248,9.869,12.32,13.84,14.77,15.34,15.89,16.06,16.05,15.94,15.78,15.58,15.37,14.24,13.17,12.2,11.36,9.946,8.829,7.929,7.191,6.576,6.057,5.613};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_12[i]*dens[j]*factor);}
G4double F_13[31]={4.243,4.737,5.669,6.544,7.369,8.148,8.886,9.587,1.25,12.93,14.66,15.76,16.46,17.18,17.45,17.51,17.44,17.31,17.13,16.93,15.81,14.71,13.71,12.81,11.31,1.09,9.109,8.292,7.606,7.024,6.524};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_13[i]*dens[j]*factor);}
G4double F_14[31]={4.384,4.889,5.844,6.742,7.591,8.397,9.164,9.895,1.59,13.47,15.41,16.68,17.5,18.4,18.78,18.9,18.88,18.78,18.63,18.45,17.35,16.22,15.19,14.25,12.66,11.36,1.3,9.409,8.658,8.016,7.462};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_14[i]*dens[j]*factor);}
G4double F_15[31]={4.538,5.054,6.03,6.951,7.825,8.656,9.45,1.21,1.93,14.01,16.15,17.59,18.56,19.65,2.16,2.36,2.41,2.35,2.23,2.07,19.03,17.9,16.84,15.86,14.18,12.79,11.64,1.67,9.844,9.133,8.517};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_15[i]*dens[j]*factor);}
G4double F_16[31]={4.67,5.2,6.202,7.148,8.047,8.904,9.723,1.51,11.26,14.52,16.86,18.47,19.58,2.88,21.52,21.82,21.93,21.92,21.84,21.71,2.73,19.61,18.53,17.52,15.76,14.29,13.05,12,11.1,1.32,9.639};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_16[i]*dens[j]*factor);}
G4double F_17[31]={4.832,5.372,6.391,7.353,8.271,9.148,9.989,1.8,11.57,14.98,17.51,19.29,2.53,22.03,22.8,23.19,23.37,23.41,23.36,23.26,22.35,21.24,2.14,19.11,17.28,15.73,14.42,13.3,12.33,11.49,1.75};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_17[i]*dens[j]*factor);}
G4double F_18[31]={4.931,5.481,6.515,7.491,8.423,9.316,1.17,11,11.79,15.33,18.02,19.96,21.33,23.03,23.93,24.41,24.64,24.73,24.71,24.64,23.79,22.68,21.58,2.52,18.64,17.03,15.66,14.49,13.47,12.58,11.8};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],F_18[i]*dens[j]*factor);}
j++;

G4double G_3[31]={1.514,1.638,1.842,2.004,2.133,2.236,2.317,2.379,2.427,2.52,2.5,2.443,2.372,2.226,2.088,1.963,1.851,1.75,1.659,1.577,1.265,1.057,.9105,.801,.6489,.5484,.4768,.4231,.3813,.3477,.32};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_3[i]*dens[j]*factor);}
G4double G_4[31]={1.913,2.085,2.375,2.615,2.815,2.982,3.122,3.237,3.331,3.582,3.636,3.612,3.556,3.409,3.253,3.102,2.96,2.829,2.707,2.594,2.143,1.824,1.587,1.406,1.148,.9733,.8473,.7521,.6777,.6179,.5687};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_4[i]*dens[j]*factor);}
G4double G_5[31]={2.238,2.453,2.825,3.138,3.406,3.637,3.836,4.007,4.152,4.594,4.753,4.786,4.762,4.642,4.488,4.327,4.169,4.018,3.874,3.738,3.169,2.743,2.415,2.158,1.779,1.516,1.323,1.176,1.06,.9664,.8893};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_5[i]*dens[j]*factor);}
G4double G_6[31]={2.506,2.761,3.209,3.592,3.925,4.217,4.474,4.7,4.898,5.549,5.838,5.948,5.97,5.897,5.761,5.603,5.44,5.278,5.12,4.968,4.304,3.783,3.369,3.034,2.529,2.169,1.9,1.692,1.527,1.394,1.283};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_6[i]*dens[j]*factor);}
G4double G_7[31]={2.732,3.021,3.534,3.981,4.375,4.724,5.036,5.315,5.563,6.432,6.87,7.074,7.158,7.154,7.051,6.909,6.75,6.586,6.422,6.261,5.528,4.925,4.431,4.022,3.39,2.927,2.576,2.301,2.08,1.9,1.75};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_7[i]*dens[j]*factor);}
G4double G_8[31]={2.934,3.251,3.824,4.329,4.779,5.183,5.547,5.876,6.173,7.259,7.855,8.168,8.324,8.406,8.348,8.231,8.086,7.927,7.764,7.599,6.818,6.147,5.582,5.105,4.349,3.783,3.345,2.998,2.717,2.485,2.291};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_8[i]*dens[j]*factor);}
G4double G_9[31]={3.114,3.456,4.078,4.635,5.135,5.589,6.001,6.376,6.718,8.011,8.77,9.196,9.431,9.608,9.602,9.514,9.386,9.236,9.076,8.911,8.097,7.374,6.75,6.214,5.349,4.687,4.168,3.751,3.41,3.127,2.888};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_9[i]*dens[j]*factor);}
G4double G_10[31]={3.28,3.643,4.309,4.911,5.457,5.956,6.412,6.831,7.216,8.707,9.63,1.18,1.5,1.78,1.84,1.79,1.68,1.54,1.39,1.23,9.397,8.63,7.956,7.367,6.401,5.649,5.05,4.565,4.163,3.827,3.542};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_10[i]*dens[j]*factor);}
G4double G_11[31]={3.422,3.803,4.507,5.148,5.736,6.277,6.775,7.236,7.661,9.347,1.44,11.12,11.54,11.95,12.08,12.07,11.99,11.87,11.73,11.58,1.75,9.946,9.227,8.591,7.53,6.689,6.011,5.454,4.99,4.598,4.264};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_11[i]*dens[j]*factor);}
G4double G_12[31]={3.568,3.963,4.697,5.371,5.994,6.572,7.107,7.604,8.067,9.931,11.19,12,12.52,13.06,13.27,13.31,13.27,13.17,13.04,12.9,12.07,11.24,1.49,9.811,8.666,7.745,6.994,6.372,5.849,5.405,5.023};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_12[i]*dens[j]*factor);}
G4double G_13[31]={3.704,4.112,4.873,5.576,6.23,6.84,7.41,7.942,8.439,1.47,11.89,12.83,13.46,14.15,14.44,14.54,14.54,14.47,14.36,14.23,13.42,12.57,11.79,11.08,9.853,8.856,8.034,7.347,6.765,6.268,5.837};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_13[i]*dens[j]*factor);}
G4double G_14[31]={3.835,4.255,5.039,5.767,6.447,7.086,7.686,8.25,8.779,1.97,12.54,13.62,14.35,15.18,15.57,15.73,15.77,15.73,15.64,15.53,14.75,13.89,13.08,12.34,11.05,9.98,9.092,8.344,7.707,7.158,6.681};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_14[i]*dens[j]*factor);}
G4double G_15[31]={3.959,4.391,5.195,5.945,6.65,7.315,7.943,8.535,9.094,11.44,13.15,14.37,15.21,16.21,16.71,16.95,17.04,17.04,16.98,16.89,16.16,15.32,14.49,13.72,12.36,11.22,1.27,9.452,8.753,8.148,7.62};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_15[i]*dens[j]*factor);}
G4double G_16[31]={4.077,4.521,5.348,6.12,6.848,7.538,8.193,8.812,9.398,11.88,13.74,15.09,16.05,17.23,17.83,18.15,18.3,18.34,18.32,18.25,17.59,16.76,15.92,15.13,13.71,12.51,11.49,1.61,9.853,9.192,8.613};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_16[i]*dens[j]*factor);}
G4double G_17[31]={4.203,4.656,5.5,6.289,7.036,7.747,8.424,9.068,9.679,12.3,14.29,15.77,16.84,18.19,18.91,19.31,19.52,19.6,19.61,19.57,18.98,18.17,17.33,16.52,15.05,13.79,12.71,11.77,1.96,1.25,9.617};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_17[i]*dens[j]*factor);}
G4double G_18[31]={4.313,4.776,5.638,6.443,7.208,7.938,8.635,9.299,9.933,12.67,14.79,16.39,17.58,19.09,19.93,2.41,2.67,2.81,2.85,2.83,2.31,19.51,18.66,17.83,16.31,15,13.86,12.88,12.02,11.26,1.59};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],G_18[i]*dens[j]*factor);}
j++;

G4double H_3[31]={1.909,2.083,2.374,2.602,2.778,2.91,3.008,3.079,3.128,3.196,3.137,3.042,2.938,2.733,2.547,2.38,2.232,2.101,1.984,1.878,1.484,1.227,1.048,.9166,.7364,.6188,.5359,.4741,.4262,.3879,.3565};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_3[i]*dens[j]*factor);}
G4double H_4[31]={2.35,2.582,2.986,3.325,3.605,3.832,4.014,4.157,4.267,4.525,4.547,4.487,4.395,4.182,3.967,3.765,3.577,3.404,3.246,3.1,2.525,2.127,1.837,1.617,1.308,1.101,.9543,.8441,.7585,.69,.6338};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_4[i]*dens[j]*factor);}
G4double H_5[31]={2.696,2.981,3.484,3.919,4.295,4.615,4.885,5.108,5.29,5.791,5.935,5.937,5.878,5.689,5.472,5.253,5.042,4.842,4.653,4.477,3.747,3.212,2.807,2.492,2.034,1.721,1.494,1.322,1.188,1.081,.9923};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_5[i]*dens[j]*factor);}
G4double H_6[31]={2.972,3.304,3.899,4.422,4.885,5.292,5.645,5.949,6.207,6.988,7.286,7.374,7.365,7.224,7.022,6.802,6.58,6.364,6.155,5.956,5.101,4.443,3.927,3.515,2.901,2.469,2.151,1.907,1.715,1.56,1.433};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_6[i]*dens[j]*factor);}
G4double H_7[31]={3.208,3.579,4.25,4.848,5.386,5.869,6.3,6.679,7.011,8.091,8.573,8.771,8.831,8.765,8.596,8.389,8.169,7.947,7.727,7.514,6.562,5.795,5.177,4.671,3.898,3.34,2.922,2.598,2.34,2.13,1.957};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_7[i]*dens[j]*factor);}
G4double H_8[31]={3.419,3.823,4.562,5.228,5.833,6.384,6.885,7.336,7.737,9.121,9.806,1.13,1.27,1.3,1.18,10,9.794,9.574,9.351,9.13,8.108,7.251,6.54,5.945,5.016,4.33,3.805,3.393,3.062,2.791,2.566};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_8[i]*dens[j]*factor);}
G4double H_9[31]={3.605,4.037,4.835,5.56,6.223,6.833,7.395,7.907,8.372,1.05,1.94,11.4,11.63,11.76,11.69,11.54,11.34,11.13,1.9,1.68,9.605,8.676,7.889,7.221,6.158,5.357,4.735,4.241,3.84,3.509,3.231};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_9[i]*dens[j]*factor);}
G4double H_10[31]={3.767,4.221,5.068,5.843,6.556,7.218,7.831,8.399,8.92,1.87,11.98,12.59,12.92,13.17,13.16,13.05,12.87,12.67,12.45,12.23,11.12,1.14,9.283,8.549,7.362,6.451,5.733,5.157,4.685,4.292,3.961};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_10[i]*dens[j]*factor);}
G4double H_11[31]={3.897,4.371,5.263,6.088,6.854,7.567,8.235,8.858,9.438,11.69,13.05,13.84,14.31,14.73,14.83,14.79,14.66,14.5,14.3,14.1,13.01,11.99,11.08,1.28,8.943,7.892,7.049,6.36,5.79,5.311,4.904};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_11[i]*dens[j]*factor);}
G4double H_12[31]={4.043,4.528,5.444,6.3,7.098,7.845,8.547,9.207,9.826,12.29,13.86,14.8,15.37,15.92,16.08,16.06,15.95,15.79,15.59,15.38,14.26,13.19,12.23,11.38,9.973,8.854,7.952,7.211,6.593,6.072,5.627};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_12[i]*dens[j]*factor);}
G4double H_13[31]={4.179,4.676,5.618,6.503,7.334,8.116,8.855,9.552,1.21,12.9,14.68,15.8,16.5,17.21,17.48,17.52,17.45,17.32,17.15,16.95,15.83,14.74,13.74,12.84,11.34,1.12,9.136,8.316,7.628,7.044,6.541};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_13[i]*dens[j]*factor);}
G4double H_14[31]={4.312,4.82,5.785,6.695,7.555,8.368,9.137,9.867,1.56,13.45,15.43,16.72,17.56,18.45,18.81,18.93,18.9,18.8,18.65,18.46,17.37,16.25,15.22,14.29,12.69,11.4,1.33,9.439,8.685,8.04,7.483};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_14[i]*dens[j]*factor);}
G4double H_15[31]={4.45,4.968,5.953,6.887,7.774,8.615,9.415,1.18,1.9,13.98,16.16,17.64,18.61,19.7,2.19,2.39,2.43,2.37,2.25,2.09,19.05,17.94,16.88,15.91,14.23,12.84,11.68,1.71,9.878,9.164,8.545};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_15[i]*dens[j]*factor);}
G4double H_16[31]={4.574,5.105,6.114,7.072,7.987,8.857,9.685,1.47,11.23,14.48,16.86,18.51,19.64,2.93,21.55,21.84,21.95,21.94,21.86,21.73,2.76,19.65,18.57,17.57,15.81,14.34,13.1,12.04,11.14,1.36,9.674};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_16[i]*dens[j]*factor);}
G4double H_17[31]={4.722,5.263,6.289,7.266,8.201,9.094,9.946,1.76,11.54,14.94,17.5,19.33,2.59,22.09,22.85,23.23,23.4,23.44,23.39,23.29,22.39,21.29,2.2,19.17,17.34,15.79,14.48,13.35,12.38,11.54,1.8};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_17[i]*dens[j]*factor);}
G4double H_18[31]={4.82,5.37,6.411,7.401,8.352,9.262,1.13,1.97,11.77,15.3,18.02,20,21.4,23.1,23.99,24.46,24.68,24.77,24.75,24.67,23.83,22.74,21.64,2.59,18.71,17.1,15.73,14.55,13.52,12.63,11.84};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],H_18[i]*dens[j]*factor);}
j++;

G4double I_3[31]={1.707,1.866,2.117,2.301,2.438,2.539,2.615,2.67,2.71,2.774,2.74,2.674,2.597,2.437,2.284,2.145,2.02,1.907,1.805,1.714,1.366,1.136,.9741,.8543,.6893,.5809,.5042,.4469,.4023,.3666,.3372};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_3[i]*dens[j]*factor);}
G4double I_4[31]={2.09,2.31,2.678,2.965,3.192,3.371,3.514,3.626,3.715,3.934,3.972,3.941,3.881,3.724,3.555,3.389,3.232,3.085,2.949,2.823,2.319,1.964,1.703,1.504,1.222,1.033,.8972,.7952,.7156,.6518,.5994};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_4[i]*dens[j]*factor);}
G4double I_5[31]={2.389,2.659,3.13,3.515,3.829,4.088,4.302,4.479,4.625,5.041,5.183,5.211,5.185,5.061,4.899,4.725,4.551,4.383,4.223,4.071,3.434,2.96,2.597,2.313,1.898,1.612,1.403,1.245,1.12,1.02,.9381};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_5[i]*dens[j]*factor);}
G4double I_6[31]={2.629,2.94,3.5,3.976,4.376,4.714,5,5.243,5.449,6.09,6.361,6.466,6.489,6.42,6.281,6.113,5.936,5.757,5.581,5.412,4.669,4.088,3.629,3.259,2.704,2.311,2.019,1.794,1.617,1.473,1.354};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_6[i]*dens[j]*factor);}
G4double I_7[31]={2.835,3.177,3.806,4.359,4.837,5.249,5.604,5.911,6.178,7.059,7.483,7.684,7.771,7.78,7.681,7.533,7.363,7.184,7.002,6.822,6.001,5.327,4.779,4.327,3.63,3.124,2.741,2.443,2.204,2.01,1.849};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_7[i]*dens[j]*factor);}
G4double I_8[31]={3.025,3.392,4.077,4.696,5.245,5.726,6.147,6.517,6.843,7.97,8.559,8.87,9.032,9.134,9.089,8.973,8.82,8.648,8.467,8.284,7.408,6.658,6.029,5.5,4.666,4.044,3.566,3.188,2.883,2.633,2.424};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_8[i]*dens[j]*factor);}
G4double I_9[31]={3.194,3.582,4.312,4.983,5.591,6.134,6.615,7.043,7.425,8.789,9.548,9.974,1.22,1.42,1.43,1.35,1.22,1.06,9.882,9.698,8.787,7.977,7.283,6.689,5.736,5.01,4.443,3.989,3.619,3.313,3.055};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_9[i]*dens[j]*factor);}
G4double I_10[31]={3.352,3.759,4.526,5.24,5.899,6.496,7.032,7.515,7.948,9.545,1.48,11.03,11.36,11.67,11.75,11.71,11.61,11.46,11.3,11.12,1.18,9.324,8.573,7.922,6.858,6.034,5.38,4.852,4.417,4.054,3.746};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_10[i]*dens[j]*factor);}
G4double I_11[31]={3.483,3.906,4.707,5.459,6.163,6.813,7.405,7.942,8.43,1.28,11.41,12.12,12.56,13.03,13.21,13.24,13.18,13.07,12.93,12.77,11.86,1.98,1.18,9.467,8.279,7.335,6.573,5.948,5.428,4.989,4.615};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_11[i]*dens[j]*factor);}
G4double I_12[31]={3.626,4.061,4.882,5.657,6.391,7.076,7.708,8.287,8.816,1.86,12.16,12.98,13.52,14.1,14.35,14.42,14.38,14.29,14.15,13.99,13.06,12.14,11.3,1.55,9.285,8.277,7.457,6.779,6.212,5.732,5.319};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_12[i]*dens[j]*factor);}
G4double I_13[31]={3.759,4.208,5.051,5.846,6.605,7.323,7.994,8.614,9.185,11.43,12.91,13.87,14.51,15.25,15.59,15.73,15.74,15.68,15.56,15.42,14.51,13.56,12.69,11.9,1.55,9.462,8.565,7.817,7.186,6.647,6.183};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_13[i]*dens[j]*factor);}
G4double I_14[31]={3.889,4.35,5.213,6.027,6.805,7.55,8.253,8.908,9.516,11.94,13.59,14.7,15.45,16.34,16.78,16.99,17.05,17.02,16.94,16.81,15.93,14.97,14.07,13.24,11.82,1.66,9.687,8.874,8.183,7.589,7.075};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_14[i]*dens[j]*factor);}
G4double I_15[31]={4.017,4.488,5.371,6.202,6.998,7.765,8.496,9.184,9.827,12.43,14.25,15.51,16.38,17.44,18,18.29,18.41,18.43,18.38,18.28,17.47,16.51,15.59,14.73,13.24,11.99,1.95,1.06,9.303,8.647,8.076};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_15[i]*dens[j]*factor);}
G4double I_16[31]={4.134,4.62,5.527,6.377,7.191,7.978,8.733,9.45,1.13,12.9,14.89,16.29,17.28,18.51,19.2,19.57,19.77,19.84,19.83,19.76,19.02,18.08,17.14,16.26,14.7,13.38,12.26,11.31,1.48,9.766,9.137};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_16[i]*dens[j]*factor);}
G4double I_17[31]={4.266,4.762,5.689,6.554,7.382,8.184,8.958,9.698,1.4,13.33,15.47,17.01,18.12,19.54,2.34,2.81,21.07,21.19,21.21,21.18,2.52,19.59,18.65,17.75,16.13,14.75,13.56,12.54,11.66,1.88,1.2};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_17[i]*dens[j]*factor);}
G4double I_18[31]={4.361,4.869,5.813,6.692,7.531,8.343,9.13,9.889,1.61,13.68,15.96,17.63,18.86,2.44,21.36,21.91,22.23,22.4,22.47,22.46,21.86,2.95,20,19.08,17.42,15.98,14.75,13.68,12.74,11.92,11.2};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],I_18[i]*dens[j]*factor);}
j++;

G4double J_3[31]={2.149,2.349,2.681,2.933,3.118,3.25,3.342,3.405,3.445,3.476,3.393,3.281,3.163,2.935,2.729,2.547,2.386,2.243,2.116,2.002,1.576,1.3,1.109,.9689,.777,.6521,.5641,.4987,.448,.4075,.3743};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_3[i]*dens[j]*factor);}
G4double J_4[31]={2.637,2.902,3.364,3.746,4.054,4.294,4.477,4.616,4.718,4.931,4.923,4.841,4.733,4.492,4.254,4.031,3.825,3.637,3.464,3.306,2.684,2.255,1.944,1.709,1.379,1.16,1.004,.8878,.7972,.7248,.6655};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_4[i]*dens[j]*factor);}
G4double J_5[31]={3.024,3.346,3.92,4.415,4.836,5.185,5.468,5.695,5.874,6.329,6.435,6.411,6.333,6.113,5.87,5.627,5.395,5.175,4.969,4.777,3.984,3.406,2.971,2.633,2.145,1.812,1.571,1.39,1.248,1.134,1.041};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_5[i]*dens[j]*factor);}
G4double J_6[31]={3.334,3.707,4.38,4.974,5.497,5.949,6.332,6.651,6.914,7.658,7.912,7.969,7.939,7.764,7.533,7.287,7.041,6.803,6.574,6.356,5.424,4.711,4.156,3.713,3.058,2.599,2.261,2.004,1.801,1.638,1.503};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_6[i]*dens[j]*factor);}
G4double J_7[31]={3.611,4.025,4.777,5.451,6.058,6.599,7.073,7.482,7.829,8.893,9.326,9.489,9.524,9.421,9.222,8.987,8.741,8.495,8.252,8.017,6.976,6.144,5.476,4.932,4.106,3.513,3.07,2.727,2.455,2.234,2.051};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_7[i]*dens[j]*factor);}
G4double J_8[31]={3.86,4.31,5.135,5.88,6.56,7.179,7.735,8.228,8.658,1.05,1.68,1.97,11.09,11.08,1.92,1.71,1.48,1.23,9.99,9.746,8.623,7.689,6.919,6.278,5.283,4.552,3.996,3.56,3.211,2.925,2.688};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_8[i]*dens[j]*factor);}
G4double J_9[31]={4.066,4.547,5.433,6.238,6.978,7.662,8.287,8.853,9.357,11.07,11.92,12.34,12.54,12.62,12.52,12.34,12.12,11.88,11.63,11.38,1.2,9.196,8.344,7.624,6.485,5.631,4.971,4.448,4.025,3.676,3.384};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_9[i]*dens[j]*factor);}
G4double J_10[31]={4.255,4.762,5.702,6.56,7.353,8.091,8.776,9.406,9.978,12.01,13.08,13.64,13.94,14.14,14.1,13.95,13.75,13.52,13.27,13.02,11.8,1.72,9.805,9.014,7.741,6.771,6.011,5.402,4.904,4.491,4.143};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_10[i]*dens[j]*factor);}
G4double J_11[31]={4.408,4.939,5.931,6.845,7.693,8.488,9.234,9.93,1.57,12.96,14.31,15.06,15.5,15.88,15.95,15.88,15.74,15.55,15.33,15.1,13.9,12.77,11.78,1.9,9.461,8.331,7.427,6.693,6.086,5.578,5.147};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_11[i]*dens[j]*factor);}
G4double J_12[31]={4.579,5.119,6.133,7.075,7.954,8.779,9.557,1.29,1.97,13.62,15.17,16.07,16.6,17.09,17.22,17.17,17.03,16.84,16.62,16.39,15.14,13.96,12.92,12,1.49,9.297,8.338,7.554,6.902,6.353,5.884};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_12[i]*dens[j]*factor);}
G4double J_13[31]={4.741,5.295,6.339,7.313,8.224,9.083,9.897,1.66,11.39,14.3,16.1,17.17,17.83,18.49,18.71,18.73,18.63,18.47,18.27,18.05,16.81,15.6,14.52,13.55,11.92,1.63,9.58,8.712,7.985,7.368,6.839};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_13[i]*dens[j]*factor);}
G4double J_14[31]={4.901,5.468,6.538,7.539,8.48,9.369,1.21,11.01,11.77,14.91,16.94,18.19,18.98,19.81,2.14,2.22,2.17,2.04,19.87,19.66,18.44,17.21,16.08,15.07,13.35,11.96,1.83,9.885,9.087,8.408,7.822};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_14[i]*dens[j]*factor);}
G4double J_15[31]={5.077,5.658,6.753,7.782,8.752,9.669,1.54,11.37,12.16,15.52,17.77,19.21,2.14,21.17,21.62,21.79,21.81,21.72,21.58,21.39,2.23,18.99,17.83,16.78,14.96,13.47,12.24,11.21,1.33,9.58,8.927};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_15[i]*dens[j]*factor);}
G4double J_16[31]={5.227,5.824,6.95,8.008,9.008,9.954,1.85,11.71,12.54,16.09,18.57,2.2,21.28,22.51,23.09,23.35,23.43,23.4,23.29,23.14,22.04,2.81,19.63,18.54,16.63,15.05,13.72,12.6,11.65,1.82,1.1};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_16[i]*dens[j]*factor);}
G4double J_17[31]={5.413,6.021,7.168,8.248,9.27,1.24,11.16,12.04,12.89,16.61,19.3,21.12,22.34,23.77,24.48,24.83,24.98,24.99,24.92,24.79,23.76,22.54,21.34,2.22,18.23,16.57,15.16,13.97,12.94,12.05,11.27};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_17[i]*dens[j]*factor);}
G4double J_18[31]={5.524,6.145,7.31,8.406,9.446,1.43,11.37,12.27,13.14,16.99,19.88,21.86,23.23,24.85,25.69,26.13,26.33,26.39,26.35,26.25,25.28,24.07,22.85,21.7,19.66,17.93,16.47,15.21,14.13,13.18,12.35};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],J_18[i]*dens[j]*factor);}
j++;

G4double K_3[31]={1.622,1.757,1.984,2.165,2.312,2.429,2.522,2.595,2.652,2.774,2.762,2.702,2.624,2.457,2.297,2.152,2.022,1.905,1.8,1.705,1.349,1.116,.9538,.8343,.6707,.5639,.4886,.4325,.3889,.354,.3255};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_3[i]*dens[j]*factor);}
G4double K_4[31]={2.038,2.224,2.541,2.806,3.029,3.218,3.376,3.507,3.616,3.918,3.997,3.981,3.922,3.757,3.577,3.402,3.237,3.084,2.943,2.813,2.294,1.933,1.67,1.471,1.19,1.003,.8695,.7696,.6918,.6295,.5785};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_4[i]*dens[j]*factor);}
G4double K_5[31]={2.368,2.601,3.005,3.348,3.645,3.903,4.127,4.32,4.486,5.003,5.207,5.26,5.241,5.109,4.932,4.746,4.561,4.384,4.217,4.059,3.403,2.918,2.551,2.264,1.85,1.566,1.36,1.205,1.083,.9854,.9052};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_5[i]*dens[j]*factor);}
G4double K_6[31]={2.629,2.905,3.391,3.81,4.178,4.503,4.791,5.046,5.27,6.024,6.377,6.522,6.559,6.483,6.326,6.142,5.95,5.76,5.575,5.397,4.628,4.032,3.565,3.192,2.636,2.245,1.957,1.736,1.562,1.422,1.306};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_6[i]*dens[j]*factor);}
G4double K_7[31]={2.856,3.166,3.721,4.207,4.639,5.025,5.373,5.686,5.966,6.963,7.487,7.744,7.854,7.86,7.74,7.571,7.382,7.187,6.993,6.803,5.947,5.254,4.695,4.237,3.538,3.034,2.655,2.362,2.129,1.939,1.782};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_7[i]*dens[j]*factor);}
G4double K_8[31]={3.056,3.396,4.013,4.561,5.053,5.498,5.902,6.269,6.603,7.841,8.547,8.93,9.126,9.234,9.166,9.024,8.849,8.658,8.462,8.265,7.346,6.572,5.928,5.39,4.55,3.929,3.455,3.083,2.784,2.539,2.335};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_8[i]*dens[j]*factor);}
G4double K_9[31]={3.218,3.582,4.251,4.853,5.398,5.894,6.349,6.765,7.148,8.61,9.495,1.01,1.3,1.51,1.5,1.39,1.23,1.05,9.855,9.656,8.694,7.857,7.147,6.544,5.584,4.861,4.299,3.853,3.491,3.191,2.941};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_9[i]*dens[j]*factor);}
G4double K_10[31]={3.364,3.749,4.461,5.11,5.702,6.246,6.747,7.21,7.637,9.314,1.38,11.03,11.42,11.76,11.81,11.74,11.6,11.43,11.24,11.04,1.05,9.166,8.398,7.737,6.667,5.846,5.199,4.68,4.254,3.9,3.602};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_10[i]*dens[j]*factor);}
G4double K_11[31]={3.476,3.878,4.632,5.327,5.97,6.565,7.116,7.629,8.107,1.03,11.31,12.14,12.66,13.18,13.35,13.35,13.27,13.13,12.97,12.79,11.83,1.9,1.08,9.351,8.14,7.185,6.419,5.793,5.275,4.84,4.47};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_11[i]*dens[j]*factor);}
G4double K_12[31]={3.604,4.015,4.786,5.505,6.174,6.798,7.38,7.924,8.434,1.52,11.96,12.91,13.53,14.18,14.41,14.44,14.37,14.24,14.08,13.9,12.89,11.93,11.07,1.31,9.035,8.027,7.213,6.544,5.987,5.517,5.115};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_12[i]*dens[j]*factor);}
G4double K_13[31]={3.723,4.144,4.937,5.679,6.378,7.034,7.649,8.228,8.772,11.03,12.64,13.75,14.49,15.31,15.65,15.75,15.72,15.62,15.48,15.31,14.32,13.33,12.43,11.63,1.27,9.178,8.286,7.547,6.926,6.398,5.945};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_13[i]*dens[j]*factor);}
G4double K_14[31]={3.84,4.27,5.081,5.844,6.566,7.248,7.893,8.502,9.077,11.5,13.27,14.52,15.39,16.38,16.83,17,17.02,16.95,16.83,16.67,15.7,14.7,13.77,12.93,11.5,1.33,9.367,8.563,7.883,7.301,6.798};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_14[i]*dens[j]*factor);}
G4double K_15[31]={3.971,4.409,5.238,6.02,6.764,7.472,8.144,8.782,9.387,11.96,13.89,15.29,16.29,17.48,18.05,18.31,18.39,18.36,18.27,18.14,17.22,16.22,15.27,14.39,12.88,11.63,1.59,9.708,8.96,8.316,7.758};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_15[i]*dens[j]*factor);}
G4double K_16[31]={4.08,4.53,5.379,6.182,6.95,7.683,8.382,9.047,9.68,12.41,14.49,16.04,17.16,18.55,19.26,19.61,19.76,19.78,19.73,19.62,18.77,17.77,16.8,15.9,14.31,12.98,11.87,1.91,1.1,9.394,8.779};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_16[i]*dens[j]*factor);}
G4double K_17[31]={4.22,4.677,5.538,6.353,7.136,7.889,8.609,9.298,9.956,12.81,15.04,16.72,17.97,19.55,2.4,2.84,21.05,21.12,21.1,21.02,2.23,19.24,18.26,17.34,15.69,14.29,13.11,12.1,11.22,1.46,9.792};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_17[i]*dens[j]*factor);}
G4double K_18[31]={4.299,4.762,5.633,6.458,7.252,8.018,8.755,9.463,1.14,13.11,15.46,17.27,18.64,2.41,21.38,21.91,22.18,22.3,22.31,22.26,21.52,2.55,19.55,18.61,16.92,15.47,14.24,13.18,12.25,11.45,1.74};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],K_18[i]*dens[j]*factor);}
j++;

G4double L_3[31]={1.968,2.15,2.444,2.663,2.824,2.943,3.029,3.091,3.133,3.185,3.123,3.027,2.922,2.717,2.53,2.363,2.215,2.083,1.966,1.861,1.466,1.211,1.033,.9028,.7245,.6083,.5265,.4656,.4184,.3807,.3498};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_3[i]*dens[j]*factor);}
G4double L_4[31]={2.403,2.652,3.075,3.412,3.679,3.89,4.055,4.184,4.283,4.511,4.526,4.464,4.371,4.157,3.942,3.739,3.55,3.377,3.218,3.072,2.497,2.1,1.811,1.593,1.286,1.083,.9375,.8289,.7446,.6771,.6219};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_4[i]*dens[j]*factor);}
G4double L_5[31]={2.739,3.045,3.58,4.026,4.396,4.701,4.952,5.157,5.324,5.779,5.908,5.906,5.846,5.656,5.438,5.218,5.006,4.804,4.615,4.437,3.706,3.171,2.768,2.454,2.001,1.691,1.467,1.298,1.166,1.06,.9734};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_5[i]*dens[j]*factor);}
G4double L_6[31]={3.003,3.358,3.992,4.538,5.005,5.403,5.74,6.024,6.262,6.98,7.255,7.335,7.324,7.181,6.978,6.756,6.533,6.314,6.104,5.904,5.045,4.386,3.872,3.462,2.853,2.426,2.112,1.872,1.683,1.531,1.405};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_6[i]*dens[j]*factor);}
G4double L_7[31]={3.237,3.628,4.341,4.97,5.522,6.003,6.421,6.781,7.092,8.093,8.541,8.726,8.781,8.712,8.541,8.332,8.109,7.884,7.662,7.447,6.488,5.72,5.102,4.599,3.832,3.28,2.868,2.548,2.294,2.088,1.918};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_7[i]*dens[j]*factor);}
G4double L_8[31]={3.448,3.869,4.647,5.349,5.977,6.536,7.03,7.464,7.845,9.135,9.774,1.08,1.22,1.24,1.12,9.934,9.723,9.5,9.274,9.051,8.019,7.157,6.446,5.853,4.931,4.251,3.733,3.327,3.002,2.735,2.514};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_8[i]*dens[j]*factor);}
G4double L_9[31]={3.63,4.073,4.905,5.665,6.357,6.982,7.543,8.043,8.489,1.06,1.89,11.33,11.55,11.67,11.6,11.45,11.25,11.03,1.81,1.58,9.495,8.562,7.775,7.109,6.053,5.26,4.646,4.158,3.764,3.438,3.166};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_9[i]*dens[j]*factor);}
G4double L_10[31]={3.799,4.263,5.137,5.946,6.692,7.375,7.996,8.558,9.064,1.91,11.94,12.52,12.83,13.07,13.06,12.94,12.76,12.55,12.33,12.1,1.98,9.992,9.139,8.407,7.228,6.326,5.619,5.051,4.587,4.202,3.877};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_10[i]*dens[j]*factor);}
G4double L_11[31]={3.933,4.415,5.328,6.184,6.983,7.724,8.407,9.032,9.603,11.75,13.03,13.78,14.23,14.64,14.74,14.69,14.57,14.4,14.2,13.99,12.89,11.86,1.95,1.15,8.813,7.768,6.931,6.25,5.686,5.213,4.812};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_11[i]*dens[j]*factor);}
G4double L_12[31]={4.086,4.576,5.507,6.386,7.213,7.988,8.709,9.376,9.991,12.35,13.82,14.72,15.26,15.79,15.94,15.92,15.81,15.64,15.44,15.23,14.08,13,12.04,11.2,9.795,8.687,7.795,7.065,6.457,5.945,5.508};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_12[i]*dens[j]*factor);}
G4double L_13[31]={4.23,4.732,5.684,6.586,7.441,8.25,9.011,9.721,1.38,12.97,14.65,15.71,16.38,17.07,17.32,17.36,17.29,17.15,16.97,16.77,15.64,14.53,13.53,12.64,11.14,9.934,8.956,8.148,7.471,6.896,6.402};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_13[i]*dens[j]*factor);}
G4double L_14[31]={4.374,4.887,5.858,6.779,7.658,8.495,9.287,1.03,1.73,13.53,15.41,16.63,17.42,18.28,18.64,18.75,18.72,18.61,18.46,18.27,17.15,16.03,14.99,14.06,12.47,11.18,1.13,9.246,8.504,7.87,7.323};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_14[i]*dens[j]*factor);}
G4double L_15[31]={4.527,5.052,6.042,6.982,7.882,8.744,9.566,1.34,11.08,14.08,16.15,17.55,18.48,19.53,2.01,2.2,2.23,2.17,2.04,19.88,18.82,17.69,16.62,15.65,13.97,12.59,11.45,1.49,9.67,8.968,8.359};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_15[i]*dens[j]*factor);}
G4double L_16[31]={4.659,5.201,6.216,7.177,8.1,8.985,9.833,1.64,11.41,14.61,16.87,18.43,19.5,2.75,21.36,21.64,21.74,21.72,21.64,21.5,2.5,19.38,18.29,17.29,15.53,14.06,12.84,11.79,1.9,1.13,9.462};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_16[i]*dens[j]*factor);}
G4double L_17[31]={4.821,5.374,6.409,7.386,8.323,9.227,1.1,1.93,11.72,15.09,17.52,19.25,2.46,21.9,22.64,23.01,23.17,23.2,23.15,23.04,22.11,2.99,19.89,18.86,17.03,15.49,14.19,13.07,12.12,11.28,1.56};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_17[i]*dens[j]*factor);}
G4double L_18[31]={4.923,5.487,6.539,7.527,8.475,9.391,1.27,11.12,11.94,15.45,18.04,19.92,21.26,22.89,23.75,24.21,24.43,24.51,24.49,24.4,23.52,22.41,21.3,2.24,18.36,16.77,15.4,14.24,13.23,12.35,11.58};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],L_18[i]*dens[j]*factor);}
j++;

G4double M_3[31]={3.507,3.775,4.211,4.539,4.78,4.952,5.072,5.152,5.203,5.217,5.066,4.872,4.671,4.291,3.955,3.664,3.41,3.188,2.992,2.819,2.187,1.79,1.517,1.32,1.052,.8799,.759,.6694,.6001,.5449,.4998};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_3[i]*dens[j]*factor);}
G4double M_4[31]={4.48,4.839,5.448,5.94,6.331,6.635,6.866,7.04,7.169,7.41,7.349,7.185,6.985,6.562,6.159,5.791,5.459,5.16,4.891,4.647,3.713,3.091,2.648,2.319,1.862,1.561,1.348,1.189,1.066,.9678,.8876};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_4[i]*dens[j]*factor);}
G4double M_5[31]={5.276,5.728,6.503,7.151,7.691,8.133,8.489,8.772,8.995,9.539,9.619,9.522,9.352,8.935,8.501,8.087,7.7,7.342,7.013,6.71,5.503,4.657,4.034,3.559,2.884,2.428,2.102,1.856,1.665,1.512,1.386};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_5[i]*dens[j]*factor);}
G4double M_6[31]={5.909,6.449,7.385,8.181,8.863,9.444,9.93,1.33,1.66,11.57,11.83,11.83,11.71,11.34,1.9,1.46,1.04,9.64,9.265,8.915,7.474,6.421,5.624,5,4.094,3.468,3.012,2.666,2.394,2.175,1.995};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_6[i]*dens[j]*factor);}
G4double M_7[31]={6.512,7.125,8.193,9.111,9.914,1.61,11.22,11.74,12.17,13.47,13.96,14.09,14.05,13.74,13.33,12.89,12.44,12.02,11.61,11.23,9.59,8.349,7.384,6.615,5.473,4.667,4.07,3.612,3.249,2.956,2.713};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_7[i]*dens[j]*factor);}
G4double M_8[31]={7.02,7.704,8.902,9.939,1.85,11.67,12.39,13.02,13.56,15.28,16.02,16.3,16.36,16.16,15.79,15.37,14.92,14.48,14.05,13.64,11.84,1.43,9.31,8.399,7.017,6.023,5.275,4.695,4.232,3.855,3.542};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_8[i]*dens[j]*factor);}
G4double M_9[31]={7.39,8.136,9.451,1.59,11.61,12.53,13.35,14.08,14.72,16.86,17.86,18.3,18.47,18.39,18.07,17.67,17.23,16.79,16.35,15.92,14.01,12.47,11.22,1.19,8.599,7.436,6.549,5.853,5.292,4.832,4.448};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_9[i]*dens[j]*factor);}
G4double M_10[31]={7.713,8.512,9.933,11.17,12.28,13.29,14.2,15.03,15.76,18.32,19.61,2.23,2.51,2.56,2.32,19.95,19.53,19.09,18.64,18.2,16.2,14.54,13.18,12.04,1.26,8.933,7.909,7.098,6.439,5.894,5.436};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_10[i]*dens[j]*factor);}
G4double M_11[31]={7.985,8.827,1.34,11.67,12.87,13.95,14.95,15.87,16.7,19.7,21.3,22.14,22.57,22.8,22.65,22.34,21.95,21.52,21.08,2.64,18.56,16.8,15.32,14.07,12.09,1.59,9.417,8.48,7.713,7.076,6.537};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_11[i]*dens[j]*factor);}
G4double M_12[31]={8.334,9.2,1.77,12.16,13.42,14.58,15.64,16.62,17.53,2.92,22.84,23.9,24.48,24.91,24.87,24.62,24.26,23.86,23.42,22.98,2.86,19,17.43,16.08,13.93,12.27,1.97,9.911,9.042,8.315,7.698};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_12[i]*dens[j]*factor);}
G4double M_13[31]={8.64,9.535,11.16,12.62,13.95,15.16,16.28,17.33,18.3,22.05,24.28,25.57,26.32,26.96,27.05,26.87,26.57,26.19,25.78,25.35,23.19,21.27,19.61,18.17,15.84,14.04,12.6,11.43,1.46,9.637,8.938};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_13[i]*dens[j]*factor);}
G4double M_14[31]={8.94,9.859,11.54,13.06,14.44,15.7,16.88,17.97,19,23.07,25.6,27.12,28.03,28.89,29.11,29.02,28.76,28.42,28.03,27.61,25.45,23.47,21.73,2.22,17.75,15.8,14.24,12.96,11.89,1.99,1.21};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_14[i]*dens[j]*factor);}
G4double M_15[31]={9.314,1.26,11.99,13.56,15,16.32,17.54,18.69,19.76,24.14,26.97,28.71,29.8,3.9,31.27,31.27,31.09,3.8,3.44,3.05,27.91,25.88,24.08,22.49,19.85,17.76,16.07,14.66,13.49,12.49,11.63};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_15[i]*dens[j]*factor);}
G4double M_16[31]={9.598,1.57,12.37,14,15.49,16.87,18.15,19.34,2.47,25.14,28.27,3.26,31.54,32.89,33.42,33.54,33.43,33.19,32.88,32.51,3.42,28.36,26.49,24.84,22.05,19.81,17.98,16.46,15.17,14.07,13.12};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_16[i]*dens[j]*factor);}
G4double M_17[31]={9.983,1.97,12.8,14.46,16,17.42,18.74,19.98,21.15,26.06,29.46,31.68,33.13,34.73,35.42,35.64,35.61,35.43,35.15,34.81,32.77,3.69,28.78,27.06,24.15,21.79,19.84,18.21,16.83,15.64,14.61};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_17[i]*dens[j]*factor);}
G4double M_18[31]={1.18,11.19,13.05,14.75,16.32,17.78,19.14,2.42,21.63,26.77,3.43,32.88,34.52,36.37,37.21,37.55,37.59,37.46,37.22,36.91,34.93,32.84,3.89,29.12,26.11,23.64,21.6,19.88,18.42,17.15,16.05};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],M_18[i]*dens[j]*factor);}
j++;

G4double N_3[31]={2.654,2.875,3.232,3.499,3.697,3.84,3.942,4.013,4.059,4.086,3.975,3.828,3.676,3.388,3.134,2.912,2.719,2.548,2.398,2.264,1.77,1.455,1.238,1.08,.8645,.7247,.6264,.5534,.4968,.4516,.4146};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_3[i]*dens[j]*factor);}
G4double N_4[31]={3.29,3.594,4.103,4.51,4.833,5.088,5.286,5.438,5.552,5.786,5.758,5.641,5.493,5.179,4.878,4.602,4.351,4.124,3.918,3.731,3.007,2.517,2.165,1.901,1.532,1.287,1.114,.984,.8833,.8027,.7368};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_4[i]*dens[j]*factor);}
G4double N_5[31]={3.785,4.166,4.817,5.355,5.801,6.169,6.471,6.715,6.912,7.417,7.519,7.465,7.347,7.046,6.727,6.42,6.131,5.863,5.614,5.384,4.456,3.794,3.301,2.922,2.377,2.007,1.74,1.539,1.381,1.255,1.152};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_5[i]*dens[j]*factor);}
G4double N_6[31]={4.166,4.615,5.399,6.062,6.626,7.106,7.512,7.852,8.137,8.959,9.23,9.267,9.2,8.94,8.626,8.306,7.995,7.699,7.419,7.155,6.057,5.238,4.609,4.112,3.382,2.873,2.499,2.214,1.99,1.81,1.661};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_6[i]*dens[j]*factor);}
G4double N_7[31]={4.507,5.009,5.902,6.671,7.339,7.919,8.421,8.854,9.225,1.39,1.86,11.02,11.02,1.84,1.55,1.23,9.915,9.603,9.302,9.015,7.779,6.82,6.063,5.452,4.532,3.876,3.386,3.008,2.708,2.465,2.264};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_7[i]*dens[j]*factor);}
G4double N_8[31]={4.804,5.352,6.342,7.21,7.973,8.648,9.242,9.765,1.22,11.74,12.44,12.73,12.83,12.74,12.5,12.2,11.89,11.57,11.26,1.95,9.606,8.524,7.649,6.928,5.82,5.012,4.399,3.92,3.536,3.223,2.962};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_8[i]*dens[j]*factor);}
G4double N_9[31]={5.032,5.618,6.689,7.641,8.489,9.247,9.924,1.53,11.06,12.92,13.86,14.3,14.49,14.51,14.33,14.05,13.75,13.43,13.11,12.8,11.38,1.2,9.227,8.415,7.144,6.199,5.471,4.896,4.431,4.048,3.727};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_9[i]*dens[j]*factor);}
G4double N_10[31]={5.239,5.856,6.994,8.019,8.942,9.776,1.53,11.21,11.82,14.01,15.19,15.8,16.1,16.25,16.12,15.89,15.6,15.29,14.97,14.65,13.17,11.91,1.85,9.956,8.532,7.456,6.616,5.946,5.399,4.945,4.562};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_10[i]*dens[j]*factor);}
G4double N_11[31]={5.425,6.072,7.277,8.376,9.379,1.29,11.13,11.89,12.58,15.17,16.66,17.51,17.97,18.33,18.33,18.18,17.95,17.67,17.38,17.07,15.58,14.24,13.08,12.07,1.44,9.172,8.168,7.356,6.687,6.127,5.653};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_11[i]*dens[j]*factor);}
G4double N_12[31]={5.628,6.282,7.507,8.639,9.68,1.64,11.52,12.33,13.07,15.91,17.62,18.62,19.19,19.67,19.72,19.59,19.36,19.09,18.79,18.47,16.92,15.53,14.32,13.28,11.57,1.24,9.174,8.308,7.589,6.985,6.471};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_12[i]*dens[j]*factor);}
G4double N_13[31]={5.824,6.493,7.751,8.919,10,11.01,11.94,12.81,13.61,16.73,18.7,19.89,2.61,21.28,21.45,21.38,21.2,2.95,2.67,2.36,18.8,17.37,16.1,14.99,13.15,11.71,1.54,9.579,8.777,8.098,7.517};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_13[i]*dens[j]*factor);}
G4double N_14[31]={6.02,6.703,7.987,9.187,1.31,11.36,12.34,13.25,14.09,17.47,19.68,21.07,21.94,22.81,23.09,23.1,22.96,22.75,22.48,22.19,2.64,19.16,17.85,16.68,14.73,13.18,11.92,1.87,9.99,9.241,8.596};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_14[i]*dens[j]*factor);}
G4double N_15[31]={6.246,6.943,8.254,9.482,1.64,11.72,12.74,13.7,14.59,18.21,2.66,22.25,23.28,24.37,24.79,24.89,24.82,24.64,24.41,24.15,22.63,21.13,19.77,18.56,16.5,14.83,13.46,12.31,11.35,1.52,9.8};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_15[i]*dens[j]*factor);}
G4double N_16[31]={6.433,7.148,8.491,9.751,1.94,12.07,13.12,14.12,15.05,18.91,21.59,23.4,24.59,25.91,26.47,26.67,26.67,26.55,26.35,26.11,24.65,23.15,21.75,2.49,18.32,16.55,15.07,13.83,12.78,11.87,11.08};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_16[i]*dens[j]*factor);}
G4double N_17[31]={6.668,7.395,8.757,1.04,11.25,12.4,13.49,14.52,15.49,19.55,22.45,24.45,25.81,27.35,28.06,28.36,28.42,28.35,28.19,27.98,26.57,25.07,23.64,22.34,2.08,18.21,16.65,15.32,14.19,13.21,12.35};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_17[i]*dens[j]*factor);}
G4double N_18[31]={6.799,7.541,8.924,1.22,11.45,12.62,13.74,14.8,15.8,2.04,23.14,25.33,26.84,28.61,29.46,29.86,29.99,29.96,29.84,29.65,28.3,26.79,25.35,24.01,21.68,19.74,18.1,16.71,15.51,14.47,13.56};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],N_18[i]*dens[j]*factor);}
j++;

G4double O_3[31]={2.784,3.026,3.417,3.703,3.906,4.047,4.141,4.202,4.238,4.231,4.103,3.948,3.79,3.492,3.23,3.001,2.801,2.625,2.469,2.331,1.821,1.496,1.272,1.108,.8863,.7423,.6413,.5662,.5082,.4618,.4239};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_3[i]*dens[j]*factor);}
G4double O_4[31]={3.441,3.767,4.322,4.765,5.11,5.374,5.57,5.714,5.818,6.004,5.951,5.822,5.666,5.341,5.031,4.746,4.487,4.252,4.039,3.845,3.096,2.589,2.225,1.952,1.571,1.319,1.14,1.007,.9035,.8208,.7532};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_4[i]*dens[j]*factor);}
G4double O_5[31]={3.957,4.361,5.063,5.649,6.133,6.524,6.834,7.078,7.266,7.715,7.782,7.711,7.583,7.27,6.941,6.625,6.327,6.05,5.792,5.554,4.592,3.905,3.395,3.001,2.438,2.057,1.782,1.575,1.413,1.284,1.178};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_5[i]*dens[j]*factor);}
G4double O_6[31]={4.362,4.834,5.668,6.385,6.997,7.514,7.942,8.292,8.576,9.344,9.568,9.581,9.501,9.227,8.903,8.573,8.253,7.947,7.657,7.384,6.245,5.394,4.741,4.226,3.47,2.944,2.559,2.266,2.036,1.85,1.698};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_6[i]*dens[j]*factor);}
G4double O_7[31]={4.74,5.265,6.205,7.026,7.747,8.374,8.912,9.367,9.748,1.87,11.29,11.41,11.4,11.19,1.89,1.57,1.24,9.915,9.604,9.307,8.022,7.024,6.237,5.602,4.65,3.971,3.466,3.077,2.769,2.519,2.313};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_7[i]*dens[j]*factor);}
G4double O_8[31]={5.063,5.635,6.671,7.587,8.405,9.133,9.773,1.33,1.81,12.3,12.94,13.2,13.27,13.16,12.91,12.6,12.27,11.95,11.62,11.31,9.911,8.784,7.871,7.122,5.972,5.136,4.503,4.009,3.615,3.293,3.025};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_8[i]*dens[j]*factor);}
G4double O_9[31]={7.121,7.892,9.331,1.67,11.93,13.12,14.24,15.31,16.32,2.61,23.59,25.53,26.82,28.26,28.93,29.21,29.27,29.19,29.03,28.81,27.36,25.79,24.3,22.94,2.59,18.64,17.02,15.65,14.48,13.47,12.59};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_9[i]*dens[j]*factor);}
G4double O_10[31]={5.543,6.188,7.373,8.442,9.413,1.3,11.11,11.84,12.49,14.73,15.84,16.4,16.66,16.77,16.63,16.38,16.08,15.76,15.44,15.11,13.56,12.25,11.15,1.22,8.743,7.63,6.764,6.074,5.512,5.046,4.655};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_10[i]*dens[j]*factor);}
G4double O_11[31]={5.732,6.403,7.645,8.776,9.811,1.76,11.64,12.45,13.18,15.8,17.2,17.93,18.32,18.59,18.53,18.33,18.06,17.76,17.44,17.11,15.53,14.14,12.95,11.93,1.3,9.044,8.057,7.261,6.609,6.064,5.603};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_11[i]*dens[j]*factor);}
G4double O_12[31]={5.969,6.654,7.93,9.104,1.19,11.19,12.12,12.99,13.78,16.75,18.42,19.34,19.87,2.29,2.33,2.19,19.95,19.67,19.36,19.04,17.42,15.98,14.72,13.63,11.86,1.48,9.379,8.487,7.749,7.129,6.601};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_12[i]*dens[j]*factor);}
G4double O_13[31]={6.184,6.887,8.201,9.414,1.54,11.59,12.57,13.48,14.33,17.62,19.56,2.69,21.35,21.96,22.1,22.03,21.84,21.58,21.29,2.98,19.36,17.87,16.55,15.39,13.48,11.98,1.78,9.788,8.963,8.266,7.67};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_13[i]*dens[j]*factor);}
G4double O_14[31]={6.399,7.118,8.463,9.71,1.87,11.96,12.98,13.94,14.83,18.4,2.6,21.93,22.73,23.53,23.79,23.79,23.64,23.42,23.15,22.85,21.24,19.71,18.33,17.12,15.1,13.49,12.18,11.1,1.2,9.43,8.769};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_14[i]*dens[j]*factor);}
G4double O_15[31]={6.653,7.389,8.766,1.05,11.24,12.37,13.42,14.42,15.36,19.19,21.66,23.19,24.15,25.16,25.55,25.63,25.56,25.38,25.14,24.87,23.3,21.74,2.32,19.05,16.91,15.18,13.76,12.58,11.58,1.73,9.995};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_15[i]*dens[j]*factor);}
G4double O_16[31]={6.858,7.616,9.031,1.35,11.58,12.74,13.84,14.87,15.85,19.93,22.66,24.41,25.53,26.76,27.29,27.48,27.47,27.34,27.15,26.9,25.39,23.82,22.36,21.05,18.79,16.94,15.42,14.14,13.05,12.11,11.3};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_16[i]*dens[j]*factor);}
G4double O_17[31]={7.121,7.892,9.331,1.67,11.93,13.12,14.24,15.31,16.32,2.61,23.59,25.53,26.82,28.26,28.93,29.21,29.27,29.19,29.03,28.81,27.36,25.79,24.3,22.94,2.59,18.64,17.02,15.65,14.48,13.47,12.59};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_17[i]*dens[j]*factor);}
G4double O_18[31]={7.278,8.066,9.53,1.89,12.17,13.38,14.53,15.62,16.66,21.14,24.34,26.49,27.93,29.6,3.41,3.78,3.91,3.88,3.76,3.57,29.17,27.6,26.09,24.69,22.25,2.22,18.52,17.08,15.84,14.77,13.83};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],O_18[i]*dens[j]*factor);}
j++;

G4double P_3[31]={1.023,1.127,1.298,1.427,1.522,1.591,1.641,1.677,1.704,1.755,1.743,1.711,1.671,1.585,1.5,1.42,1.346,1.279,1.218,1.162,.9451,.7987,.6937,.6149,.5043,.4302,.3768,.3364,.3046,.2788,.2575};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_3[i]*dens[j]*factor);}
G4double P_4[31]={1.245,1.386,1.628,1.825,1.983,2.106,2.202,2.276,2.335,2.487,2.524,2.517,2.489,2.41,2.319,2.227,2.138,2.053,1.972,1.897,1.589,1.367,1.201,1.073,.8873,.7599,.667,.596,.54,.4945,.4569};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_4[i]*dens[j]*factor);}
G4double P_5[31]={1.418,1.587,1.888,2.146,2.363,2.543,2.689,2.808,2.904,3.184,3.29,3.323,3.32,3.266,3.184,3.091,2.996,2.9,2.808,2.719,2.336,2.043,1.816,1.635,1.366,1.177,1.036,.9275,.8412,.7709,.7124};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_5[i]*dens[j]*factor);}
G4double P_6[31]={1.557,1.747,2.097,2.406,2.676,2.91,3.107,3.273,3.411,3.837,4.029,4.114,4.146,4.132,4.07,3.986,3.892,3.794,3.695,3.598,3.158,2.803,2.517,2.285,1.93,1.674,1.48,1.328,1.207,1.107,1.024};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_6[i]*dens[j]*factor);}
G4double P_7[31]={1.678,1.883,2.268,2.618,2.933,3.213,3.458,3.669,3.85,4.436,4.726,4.876,4.951,4.991,4.96,4.894,4.809,4.715,4.616,4.515,4.036,3.63,3.292,3.011,2.572,2.246,1.995,1.796,1.635,1.502,1.391};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_7[i]*dens[j]*factor);}
G4double P_8[31]={1.791,2.007,2.419,2.803,3.155,3.476,3.764,4.019,4.242,4.994,5.393,5.614,5.74,5.845,5.854,5.813,5.744,5.658,5.563,5.463,4.961,4.512,4.129,3.802,3.281,2.886,2.576,2.328,2.125,1.955,1.812};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_8[i]*dens[j]*factor);}
G4double P_9[31]={1.889,2.115,2.548,2.958,3.34,3.694,4.02,4.313,4.575,5.493,6.004,6.304,6.485,6.663,6.716,6.704,6.654,6.58,6.493,6.397,5.885,5.404,4.984,4.62,4.026,3.567,3.201,2.905,2.659,2.453,2.278};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_9[i]*dens[j]*factor);}
G4double P_10[31]={1.981,2.217,2.667,3.097,3.504,3.885,4.241,4.569,4.867,5.948,6.577,6.96,7.202,7.46,7.564,7.585,7.559,7.501,7.424,7.335,6.824,6.321,5.868,5.47,4.811,4.29,3.871,3.527,3.239,2.996,2.788};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_10[i]*dens[j]*factor);}
G4double P_11[31]={2.062,2.309,2.776,3.225,3.654,4.061,4.446,4.806,5.138,6.389,7.15,7.631,7.948,8.31,8.483,8.553,8.561,8.529,8.471,8.395,7.906,7.387,6.906,6.472,5.739,5.149,4.665,4.264,3.925,3.637,3.388};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_11[i]*dens[j]*factor);}
G4double P_12[31]={2.147,2.4,2.877,3.336,3.778,4.201,4.604,4.986,5.344,6.735,7.608,8.174,8.553,9,9.226,9.334,9.368,9.354,9.308,9.241,8.762,8.231,7.73,7.275,6.498,5.865,5.342,4.905,4.532,4.213,3.937};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_12[i]*dens[j]*factor);}
G4double P_13[31]={2.227,2.489,2.98,3.448,3.903,4.341,4.761,5.163,5.544,7.074,8.069,8.729,9.181,9.73,1.02,1.18,1.25,1.26,1.24,1.19,9.739,9.205,8.686,8.209,7.381,6.698,6.128,5.645,5.233,4.877,4.566};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_13[i]*dens[j]*factor);}
G4double P_14[31]={2.303,2.574,3.079,3.557,4.021,4.47,4.905,5.322,5.722,7.378,8.491,9.246,9.773,1.43,1.79,11,11.11,11.15,11.15,11.11,1.71,1.17,9.643,9.147,8.275,7.546,6.931,6.407,5.957,5.565,5.222};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_14[i]*dens[j]*factor);}
G4double P_15[31]={2.379,2.659,3.178,3.665,4.137,4.597,5.043,5.474,5.89,7.66,8.89,9.743,1.35,11.12,11.56,11.83,11.98,12.06,12.08,12.07,11.72,11.2,1.66,1.15,9.236,8.462,7.802,7.235,6.743,6.314,5.935};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_15[i]*dens[j]*factor);}
G4double P_16[31]={1.791,2.007,2.419,2.803,3.155,3.476,3.764,4.019,4.242,4.994,5.393,5.614,5.74,5.845,5.854,5.813,5.744,5.658,5.563,5.463,4.961,4.512,4.129,3.802,3.281,2.886,2.576,2.328,2.125,1.955,1.812};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_16[i]*dens[j]*factor);}
G4double P_17[31]={2.525,2.82,3.37,3.879,4.368,4.844,5.308,5.761,6.202,8.163,9.615,1.66,11.42,12.43,13.04,13.43,13.68,13.83,13.92,13.95,13.74,13.26,12.72,12.19,11.21,1.35,9.614,8.967,8.4,7.898,7.452};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_17[i]*dens[j]*factor);}
G4double P_18[31]={2.583,2.886,3.449,3.97,4.466,4.948,5.419,5.879,6.328,8.364,9.915,11.05,11.89,13.02,13.71,14.16,14.46,14.65,14.76,14.82,14.67,14.22,13.68,13.15,12.14,11.26,1.48,9.805,9.206,8.674,8.2};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],P_18[i]*dens[j]*factor);}
j++;

G4double Q_3[31]={2.597,2.831,3.212,3.493,3.691,3.828,3.919,3.977,4.012,4.006,3.89,3.747,3.601,3.326,3.081,2.867,2.679,2.513,2.366,2.235,1.75,1.439,1.224,1.068,.8543,.7158,.6185,.5463,.4904,.4457,.4092};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_3[i]*dens[j]*factor);}
G4double Q_4[31]={3.198,3.511,4.048,4.481,4.821,5.078,5.269,5.409,5.509,5.688,5.643,5.528,5.387,5.088,4.801,4.536,4.294,4.073,3.872,3.689,2.977,2.493,2.144,1.881,1.515,1.273,1.1,.9718,.8721,.7924,.7272};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_4[i]*dens[j]*factor);}
G4double Q_5[31]={3.674,4.056,4.73,5.3,5.774,6.157,6.461,6.698,6.881,7.313,7.382,7.323,7.21,6.926,6.626,6.333,6.056,5.797,5.555,5.331,4.419,3.763,3.273,2.896,2.354,1.985,1.72,1.52,1.365,1.24,1.137};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_5[i]*dens[j]*factor);}
G4double Q_6[31]={4.051,4.495,5.289,5.98,6.577,7.083,7.502,7.845,8.121,8.862,9.081,9.103,9.036,8.794,8.5,8.198,7.901,7.617,7.346,7.09,6.012,5.2,4.574,4.079,3.351,2.844,2.472,2.189,1.967,1.787,1.64};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_6[i]*dens[j]*factor);}
G4double Q_7[31]={4.399,4.891,5.779,6.565,7.264,7.876,8.403,8.849,9.221,1.31,1.71,1.84,1.84,1.67,1.4,1.11,9.804,9.506,9.216,8.938,7.725,6.775,6.021,5.411,4.493,3.838,3.351,2.975,2.677,2.435,2.235};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_7[i]*dens[j]*factor);}
G4double Q_8[31]={4.708,5.245,6.219,7.091,7.878,8.585,9.21,9.755,1.22,11.68,12.29,12.55,12.63,12.55,12.33,12.05,11.76,11.46,11.16,1.87,9.548,8.476,7.603,6.883,5.775,4.967,4.355,3.878,3.496,3.184,2.926};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_8[i]*dens[j]*factor);}
G4double Q_9[31]={4.95,5.523,6.572,7.516,8.376,9.161,9.871,1.5,11.06,12.88,13.72,14.1,14.27,14.29,14.12,13.87,13.59,13.29,12.99,12.69,11.3,1.13,9.163,8.353,7.083,6.139,5.414,4.841,4.378,3.998,3.68};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_9[i]*dens[j]*factor);}
G4double Q_10[31]={5.171,5.776,6.888,7.895,8.819,9.67,1.45,11.16,11.8,13.98,15.06,15.6,15.86,15.99,15.88,15.67,15.4,15.11,14.81,14.51,13.06,11.82,1.76,9.872,8.452,7.378,6.542,5.875,5.331,4.881,4.502};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_10[i]*dens[j]*factor);}
G4double Q_11[31]={5.362,5.997,7.176,8.252,9.245,1.17,11.03,11.82,12.54,15.16,16.55,17.3,17.72,18.06,18.08,17.95,17.74,17.49,17.22,16.93,15.49,14.17,13.02,12.02,1.39,9.122,8.117,7.303,6.634,6.075,5.602};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_11[i]*dens[j]*factor);}
G4double Q_12[31]={5.572,6.215,7.413,8.517,9.541,1.5,11.39,12.23,13,15.9,17.51,18.4,18.91,19.35,19.41,19.3,19.1,18.85,18.57,18.27,16.77,15.4,14.21,13.17,11.46,1.13,9.076,8.213,7.499,6.899,6.388};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_12[i]*dens[j]*factor);}
G4double Q_13[31]={5.774,6.434,7.668,8.809,9.871,1.87,11.8,12.68,13.51,16.71,18.59,19.68,2.32,2.93,21.1,21.05,2.9,2.67,2.41,2.13,18.63,17.22,15.96,14.86,13.03,11.59,1.43,9.471,8.674,7.999,7.422};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_13[i]*dens[j]*factor);}
G4double Q_14[31]={5.976,6.651,7.915,9.089,1.18,11.21,12.18,13.1,13.97,17.45,19.59,2.86,21.64,22.43,22.7,22.73,22.62,22.43,22.19,21.92,2.43,18.99,17.69,16.53,14.59,13.05,11.79,1.74,9.87,9.126,8.486};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_14[i]*dens[j]*factor);}
G4double Q_15[31]={6.211,6.903,8.199,9.405,1.53,11.59,12.6,13.55,14.45,18.19,2.58,22.06,22.99,23.98,24.38,24.49,24.45,24.3,24.1,23.86,22.42,2.95,19.61,18.4,16.34,14.68,13.32,12.18,11.21,1.39,9.676};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_15[i]*dens[j]*factor);}
G4double Q_16[31]={6.402,7.115,8.448,9.69,1.85,11.95,12.98,13.97,14.91,18.88,21.54,23.22,24.31,25.51,26.05,26.26,26.28,26.19,26.02,25.8,24.43,22.96,21.58,2.33,18.17,16.39,14.92,13.69,12.63,11.73,1.94};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_16[i]*dens[j]*factor);}
G4double Q_17[31]={6.647,7.373,8.73,9.996,11.18,12.3,13.37,14.38,15.35,19.51,22.41,24.29,25.53,26.94,27.61,27.91,28,27.96,27.83,27.64,26.33,24.86,23.46,22.16,19.91,18.04,16.48,15.16,14.03,13.05,12.2};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_17[i]*dens[j]*factor);}
G4double Q_18[31]={6.781,7.522,8.903,1.19,11.4,12.54,13.62,14.65,15.64,19.98,23.1,25.17,26.56,28.17,28.98,29.37,29.52,29.52,29.42,29.26,28.01,26.54,25.12,23.79,21.47,19.53,17.9,16.51,15.32,14.28,13.38};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],Q_18[i]*dens[j]*factor);}
j++;

G4double R_3[31]={2.358,2.573,2.927,3.191,3.383,3.517,3.609,3.671,3.709,3.726,3.628,3.501,3.37,3.119,2.895,2.698,2.524,2.37,2.234,2.112,1.658,1.366,1.163,1.015,.8137,.6824,.59,.5214,.4682,.4258,.391};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_3[i]*dens[j]*factor);}
G4double R_4[31]={2.899,3.186,3.682,4.086,4.408,4.656,4.843,4.983,5.085,5.286,5.263,5.166,5.041,4.773,4.512,4.269,4.046,3.842,3.656,3.486,2.822,2.367,2.038,1.79,1.443,1.213,1.05,.9278,.8329,.7571,.695};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_4[i]*dens[j]*factor);}
G4double R_5[31]={3.326,3.677,4.296,4.825,5.269,5.633,5.926,6.158,6.34,6.789,6.88,6.841,6.747,6.497,6.226,5.959,5.706,5.467,5.245,5.037,4.188,3.574,3.113,2.757,2.243,1.894,1.642,1.452,1.303,1.184,1.087};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_5[i]*dens[j]*factor);}
G4double R_6[31]={3.666,4.073,4.802,5.44,5.996,6.471,6.871,7.201,7.471,8.219,8.46,8.503,8.456,8.249,7.988,7.715,7.445,7.185,6.936,6.7,5.699,4.94,4.352,3.886,3.196,2.715,2.361,2.092,1.88,1.709,1.568};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_6[i]*dens[j]*factor);}
G4double R_7[31]={3.974,4.425,5.241,5.967,6.615,7.188,7.685,8.11,8.47,9.55,9.975,1.12,1.14,10,9.776,9.513,9.24,8.969,8.704,8.448,7.327,6.439,5.732,5.158,4.289,3.667,3.203,2.844,2.56,2.329,2.139};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_7[i]*dens[j]*factor);}
G4double R_8[31]={4.249,4.74,5.635,6.439,7.168,7.827,8.414,8.93,9.377,1.8,11.43,11.71,11.81,11.77,11.58,11.34,11.08,1.8,1.53,1.26,9.055,8.057,7.239,6.562,5.515,4.748,4.166,3.711,3.346,3.048,2.801};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_8[i]*dens[j]*factor);}
G4double R_9[31]={4.469,4.993,5.956,6.826,7.623,8.352,9.016,9.611,1.13,11.9,12.75,13.16,13.35,13.4,13.27,13.06,12.81,12.54,12.26,11.99,1.71,9.634,8.728,7.967,6.767,5.871,5.181,4.635,4.193,3.829,3.525};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_9[i]*dens[j]*factor);}
G4double R_10[31]={4.672,5.224,6.244,7.173,8.027,8.818,9.547,1.21,1.81,12.92,14,14.55,14.84,15.01,14.93,14.75,14.52,14.26,13.99,13.72,12.39,11.23,1.25,9.418,8.077,7.058,6.263,5.627,5.108,4.677,4.314};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_10[i]*dens[j]*factor);}
G4double R_11[31]={4.841,5.42,6.5,7.49,8.406,9.261,1.06,1.8,11.47,13.97,15.35,16.11,16.54,16.9,16.95,16.85,16.68,16.46,16.21,15.96,14.64,13.43,12.36,11.42,9.896,8.702,7.752,6.981,6.345,5.813,5.362};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_11[i]*dens[j]*factor);}
G4double R_12[31]={5.029,5.617,6.716,7.735,8.682,9.568,1.4,11.18,11.9,14.66,16.25,17.16,17.68,18.15,18.25,18.17,18,17.78,17.53,17.27,15.9,14.64,13.52,12.55,1.95,9.694,8.688,7.867,7.186,6.613,6.125};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_12[i]*dens[j]*factor);}
G4double R_13[31]={5.209,5.812,6.943,7.995,8.978,9.901,1.77,11.59,12.36,15.41,17.25,18.34,18.99,19.64,19.83,19.82,19.7,19.51,19.28,19.02,17.66,16.36,15.19,14.16,12.45,11.08,9.982,9.073,8.313,7.669,7.118};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_13[i]*dens[j]*factor);}
G4double R_14[31]={5.387,6.004,7.163,8.245,9.259,1.21,11.11,11.97,12.78,16.07,18.16,19.43,2.22,21.04,21.34,21.4,21.32,21.16,2.96,2.72,19.37,18.04,16.83,15.76,13.93,12.47,11.28,1.29,9.46,8.75,8.139};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_14[i]*dens[j]*factor);}
G4double R_15[31]={5.591,6.222,7.409,8.52,9.564,1.54,11.48,12.37,13.21,16.74,19.07,2.53,21.47,22.49,22.92,23.06,23.05,22.94,22.76,22.55,21.25,19.91,18.67,17.54,15.61,14.04,12.75,11.66,1.75,9.965,9.285};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_15[i]*dens[j]*factor);}
G4double R_16[31]={5.759,6.409,7.629,8.771,9.848,1.86,11.82,12.74,13.62,17.37,19.93,21.6,22.69,23.91,24.48,24.72,24.77,24.71,24.57,24.39,23.16,21.82,2.55,19.38,17.35,15.68,14.29,13.12,12.11,11.25,1.5};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_16[i]*dens[j]*factor);}
G4double R_17[31]={5.972,6.633,7.875,9.04,1.14,11.18,12.17,13.11,14.02,17.94,2.73,22.58,23.82,25.25,25.95,26.28,26.4,26.38,26.28,26.13,24.97,23.63,22.33,21.13,19.02,17.26,15.79,14.53,13.45,12.52,11.71};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_17[i]*dens[j]*factor);}
G4double R_18[31]={6.093,6.768,8.031,9.213,1.33,11.39,12.39,13.36,14.28,18.36,21.35,23.39,24.78,26.41,27.24,27.65,27.83,27.86,27.8,27.67,26.56,25.23,23.92,22.68,2.52,18.69,17.14,15.83,14.69,13.71,12.84};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],R_18[i]*dens[j]*factor);}
j++;

G4double S_3[31]={3.038,3.3,3.719,4.021,4.232,4.375,4.469,4.527,4.56,4.531,4.384,4.211,4.037,3.712,3.428,3.181,2.966,2.777,2.61,2.462,1.918,1.574,1.337,1.164,.9296,.778,.6717,.5928,.5318,.4832,.4433};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_3[i]*dens[j]*factor);}
G4double S_4[31]={3.761,4.114,4.711,5.183,5.547,5.82,6.02,6.164,6.267,6.432,6.358,6.209,6.036,5.678,5.34,5.03,4.751,4.498,4.269,4.062,3.261,2.722,2.337,2.048,1.647,1.382,1.194,1.054,.9453,.8586,.7877};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_4[i]*dens[j]*factor);}
G4double S_5[31]={4.329,4.766,5.525,6.153,6.667,7.077,7.399,7.648,7.839,8.273,8.318,8.226,8.079,7.729,7.368,7.023,6.7,6.4,6.123,5.866,4.837,4.106,3.565,3.149,2.555,2.154,1.865,1.647,1.478,1.342,1.231};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_5[i]*dens[j]*factor);}
G4double S_6[31]={4.774,5.286,6.188,6.959,7.613,8.16,8.609,8.972,9.264,1.03,1.23,1.22,1.12,9.809,9.45,9.088,8.738,8.406,8.093,7.799,6.576,5.669,4.976,4.431,3.634,3.081,2.677,2.369,2.128,1.934,1.774};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_6[i]*dens[j]*factor);}
G4double S_7[31]={5.194,5.763,6.779,7.664,8.436,9.103,9.671,1.15,1.54,11.67,12.07,12.17,12.14,11.9,11.56,11.2,1.84,1.49,1.15,9.827,8.445,7.38,6.543,5.872,4.867,4.154,3.624,3.216,2.893,2.632,2.416};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_7[i]*dens[j]*factor);}
G4double S_8[31]={5.553,6.174,7.294,8.281,9.159,9.936,1.61,11.2,11.7,13.23,13.85,14.09,14.14,13.99,13.7,13.35,13,12.64,12.28,11.94,1.43,9.228,8.257,7.463,6.249,5.369,4.705,4.188,3.775,3.438,3.158};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_8[i]*dens[j]*factor);}
G4double S_9[31]={5.826,6.49,7.698,8.774,9.741,1.61,11.39,12.07,12.67,14.59,15.45,15.83,15.97,15.92,15.68,15.36,15.01,14.65,14.3,13.94,12.34,11.03,9.95,9.055,7.662,6.633,5.845,5.225,4.725,4.313,3.97};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_9[i]*dens[j]*factor);}
G4double S_10[31]={6.073,6.774,8.058,9.21,1.25,11.21,12.07,12.84,13.53,15.85,16.96,17.5,17.75,17.82,17.64,17.35,17.02,16.66,16.3,15.94,14.27,12.87,11.69,1.7,9.143,7.971,7.062,6.339,5.751,5.264,4.855};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_10[i]*dens[j]*factor);}
G4double S_11[31]={6.28,7.01,8.355,9.574,1.69,11.71,12.65,13.51,14.28,17.02,18.42,19.15,19.52,19.75,19.66,19.42,19.11,18.77,18.42,18.06,16.34,14.85,13.58,12.5,1.77,9.448,8.411,7.577,6.894,6.324,5.843};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_11[i]*dens[j]*factor);}
G4double S_12[31]={6.544,7.288,8.67,9.936,11.1,12.18,13.17,14.1,14.94,18.04,19.74,2.67,21.17,21.57,21.57,21.39,21.12,2.8,2.45,2.1,18.35,16.78,15.44,14.28,12.4,1.95,9.793,8.857,8.083,7.434,6.883};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_12[i]*dens[j]*factor);}
G4double S_13[31]={6.783,7.547,8.969,1.28,11.49,12.61,13.66,14.64,15.54,18.99,2.98,22.11,22.76,23.34,23.45,23.34,23.11,22.82,22.5,22.15,2.39,18.77,17.36,16.13,14.11,12.52,11.25,1.21,9.349,8.619,7.996};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_13[i]*dens[j]*factor);}
G4double S_14[31]={7.021,7.803,9.26,1.61,11.86,13.02,14.11,15.13,16.09,19.85,22.11,23.44,24.24,25.01,25.24,25.2,25.02,24.76,24.46,24.13,22.37,2.71,19.24,17.94,15.8,14.1,12.72,11.59,1.64,9.832,9.14};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_14[i]*dens[j]*factor);}
G4double S_15[31]={7.309,8.11,9.602,1.98,12.27,13.47,14.6,15.67,16.67,2.71,23.26,24.8,25.76,26.75,27.11,27.16,27.05,26.84,26.57,26.26,24.53,22.85,21.32,19.97,17.68,15.86,14.36,13.12,12.08,11.18,1.41};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_15[i]*dens[j]*factor);}
G4double S_16[31]={7.538,8.363,9.898,11.32,12.65,13.89,15.06,16.17,17.21,21.53,24.36,26.12,27.25,28.46,28.96,29.12,29.08,28.92,28.68,28.4,26.74,25.03,23.47,22.05,19.65,17.7,16.09,14.74,13.6,12.62,11.77};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_16[i]*dens[j]*factor);}
G4double S_17[31]={7.835,8.674,1.23,11.68,13.04,14.31,15.51,16.65,17.73,22.27,25.36,27.34,28.63,3.06,3.7,3.95,3.98,3.87,3.67,3.42,28.81,27.1,25.5,24.04,21.53,19.47,17.76,16.32,15.09,14.03,13.11};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_17[i]*dens[j]*factor);}
G4double S_18[31]={8.007,8.865,1.45,11.92,13.3,14.6,15.82,16.99,18.1,22.84,26.18,28.37,29.83,31.48,32.27,32.61,32.71,32.65,32.49,32.26,3.71,28.99,27.36,25.86,23.27,21.12,19.33,17.81,16.51,15.38,14.4};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],S_18[i]*dens[j]*factor);}
j++;

G4double T_3[31]={2.242,2.448,2.786,3.041,3.229,3.363,3.457,3.52,3.562,3.591,3.503,3.385,3.26,3.021,2.806,2.617,2.45,2.302,2.17,2.052,1.613,1.331,1.134,.9906,.7941,.6663,.5764,.5095,.4576,.4162,.3823};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_3[i]*dens[j]*factor);}
G4double T_4[31]={2.756,3.03,3.503,3.891,4.203,4.446,4.632,4.772,4.877,5.093,5.081,4.993,4.877,4.622,4.373,4.14,3.926,3.73,3.552,3.388,2.746,2.306,1.987,1.746,1.409,1.185,1.025,.9067,.8141,.7401,.6796};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_4[i]*dens[j]*factor);}
G4double T_5[31]={3.161,3.495,4.086,4.592,5.019,5.372,5.66,5.889,6.072,6.535,6.639,6.611,6.525,6.291,6.033,5.779,5.536,5.308,5.094,4.894,4.076,3.482,3.035,2.689,2.19,1.85,1.604,1.419,1.274,1.158,1.063};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_5[i]*dens[j]*factor);}
G4double T_6[31]={3.482,3.87,4.566,5.176,5.709,6.168,6.556,6.88,7.147,7.905,8.161,8.215,8.177,7.987,7.741,7.482,7.224,6.975,6.737,6.51,5.547,4.814,4.244,3.791,3.121,2.653,2.308,2.045,1.838,1.671,1.534};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_6[i]*dens[j]*factor);}
G4double T_7[31]={3.769,4.2,4.981,5.675,6.296,6.847,7.328,7.742,8.095,9.177,9.617,9.778,9.808,9.69,9.475,9.225,8.966,8.707,8.454,8.209,7.131,6.275,5.59,5.034,4.19,3.584,3.132,2.782,2.504,2.279,2.093};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_7[i]*dens[j]*factor);}
G4double T_8[31]={4.026,4.495,5.352,6.122,6.821,7.453,8.019,8.519,8.955,1.37,11.02,11.3,11.41,11.39,11.22,11,1.75,1.49,1.23,9.978,8.813,7.851,7.061,6.405,5.388,4.642,4.074,3.63,3.274,2.983,2.742};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_8[i]*dens[j]*factor);}
G4double T_9[31]={4.235,4.735,5.656,6.491,7.254,7.955,8.593,9.168,9.679,11.42,12.28,12.71,12.9,12.98,12.86,12.67,12.43,12.18,11.91,11.65,1.43,9.391,8.515,7.778,6.613,5.742,5.069,4.536,4.104,3.749,3.451};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_9[i]*dens[j]*factor);}
G4double T_10[31]={4.427,4.954,5.93,6.821,7.641,8.4,9.1,9.742,1.32,12.39,13.48,14.05,14.34,14.53,14.48,14.32,14.1,13.85,13.6,13.33,12.07,1.95,10,9.197,7.895,6.904,6.129,5.508,5.001,4.58,4.225};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_10[i]*dens[j]*factor);}
G4double T_11[31]={4.585,5.138,6.169,7.118,7.997,8.817,9.583,1.29,1.94,13.38,14.75,15.53,15.97,16.35,16.41,16.32,16.16,15.95,15.72,15.48,14.22,13.06,12.03,11.13,9.654,8.497,7.574,6.824,6.205,5.686,5.247};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_11[i]*dens[j]*factor);}
G4double T_12[31]={4.762,5.323,6.375,7.352,8.263,9.115,9.915,1.66,11.36,14.05,15.64,16.55,17.09,17.58,17.69,17.63,17.48,17.27,17.04,16.79,15.48,14.27,13.19,12.25,1.7,9.481,8.502,7.701,7.036,6.476,5.999};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_12[i]*dens[j]*factor);}
G4double T_13[31]={4.93,5.505,6.587,7.597,8.542,9.43,1.26,11.06,11.8,14.76,16.59,17.68,18.35,19.01,19.23,19.23,19.12,18.94,18.73,18.49,17.19,15.94,14.82,13.83,12.16,1.84,9.767,8.881,8.14,7.511,6.972};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_13[i]*dens[j]*factor);}
G4double T_14[31]={5.097,5.685,6.793,7.83,8.806,9.725,1.59,11.42,12.2,15.39,17.45,18.73,19.54,2.37,2.69,2.77,2.7,2.56,2.36,2.14,18.86,17.59,16.42,15.38,13.62,12.2,11.04,1.07,9.264,8.571,7.973};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_14[i]*dens[j]*factor);}
G4double T_15[31]={5.283,5.884,7.018,8.083,9.087,1.03,1.93,11.79,12.61,16.02,18.32,19.79,2.74,21.77,22.22,22.38,22.38,22.28,22.12,21.92,2.69,19.41,18.21,17.12,15.26,13.74,12.48,11.42,1.53,9.763,9.098};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_15[i]*dens[j]*factor);}
G4double T_16[31]={5.44,6.058,7.223,8.317,9.352,1.33,11.26,12.15,12.99,16.62,19.14,2.8,21.9,23.15,23.73,23.98,24.05,24,23.88,23.71,22.54,21.26,2.04,18.91,16.96,15.34,13.99,12.84,11.87,11.02,1.29};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_16[i]*dens[j]*factor);}
G4double T_17[31]={5.636,6.266,7.45,8.566,9.623,1.62,11.58,12.49,13.36,17.16,19.89,21.74,22.99,24.44,25.16,25.5,25.63,25.63,25.54,25.4,24.3,23.03,21.79,2.63,18.59,16.89,15.45,14.23,13.18,12.27,11.48};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_17[i]*dens[j]*factor);}
G4double T_18[31]={5.751,6.393,7.597,8.728,9.803,1.82,11.79,12.73,13.62,17.56,2.49,22.52,23.91,25.56,26.4,26.83,27.02,27.07,27.02,26.9,25.86,24.59,23.33,22.15,2.05,18.28,16.78,15.5,14.4,13.43,12.59};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],T_18[i]*dens[j]*factor);}
j++;

G4double U_3[31]={2.534,2.757,3.121,3.391,3.587,3.725,3.82,3.884,3.924,3.937,3.829,3.69,3.547,3.276,3.035,2.824,2.639,2.476,2.331,2.202,1.725,1.42,1.208,1.054,.8438,.7073,.6114,.5401,.485,.4409,.4048};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_3[i]*dens[j]*factor);}
G4double U_4[31]={3.127,3.427,3.941,4.355,4.683,4.936,5.128,5.272,5.378,5.583,5.551,5.442,5.304,5.011,4.728,4.466,4.228,4.011,3.814,3.634,2.934,2.458,2.114,1.857,1.496,1.257,1.088,.9608,.8624,.7838,.7194};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_4[i]*dens[j]*factor);}
G4double U_5[31]={3.592,3.963,4.61,5.155,5.609,5.981,6.28,6.518,6.705,7.168,7.256,7.206,7.097,6.819,6.523,6.234,5.961,5.706,5.468,5.248,4.352,3.709,3.228,2.857,2.324,1.961,1.7,1.503,1.349,1.226,1.125};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_5[i]*dens[j]*factor);}
G4double U_6[31]={3.958,4.39,5.159,5.822,6.393,6.88,7.288,7.625,7.902,8.673,8.918,8.952,8.892,8.656,8.367,8.068,7.776,7.496,7.23,6.978,5.92,5.124,4.51,4.024,3.308,2.809,2.443,2.164,1.945,1.768,1.623};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_6[i]*dens[j]*factor);}
G4double U_7[31]={4.294,4.776,5.641,6.401,7.071,7.659,8.167,8.602,8.969,1.08,1.51,1.66,1.67,1.5,1.24,9.947,9.649,9.355,9.07,8.796,7.607,6.675,5.936,5.338,4.436,3.792,3.312,2.941,2.647,2.409,2.212};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_7[i]*dens[j]*factor);}
G4double U_8[31]={4.583,5.108,6.06,6.907,7.666,8.345,8.947,9.474,9.931,11.4,12.05,12.32,12.41,12.34,12.13,11.86,11.57,11.27,1.98,1.69,9.398,8.348,7.493,6.788,5.701,4.907,4.305,3.835,3.458,3.151,2.895};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_8[i]*dens[j]*factor);}
G4double U_9[31]={4.813,5.373,6.399,7.321,8.155,8.912,9.594,1.2,1.74,12.56,13.43,13.85,14.03,14.06,13.9,13.65,13.37,13.08,12.78,12.49,11.12,9.983,9.035,8.241,6.995,6.067,5.353,4.789,4.333,3.957,3.643};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_9[i]*dens[j]*factor);}
G4double U_10[31]={5.022,5.613,6.702,7.688,8.588,9.413,1.17,1.85,11.47,13.63,14.74,15.31,15.59,15.74,15.64,15.43,15.17,14.88,14.59,14.29,12.87,11.65,1.62,9.745,8.351,7.295,6.472,5.814,5.278,4.833,4.458};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_10[i]*dens[j]*factor);}
G4double U_11[31]={5.194,5.808,6.949,7.992,8.951,9.836,1.65,11.41,12.1,14.61,15.99,16.74,17.14,17.44,17.42,17.26,17.03,16.76,16.47,16.18,14.73,13.44,12.33,11.38,9.835,8.647,7.709,6.952,6.329,5.808,5.368};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_11[i]*dens[j]*factor);}
G4double U_12[31]={5.404,6.031,7.205,8.289,9.292,1.22,11.09,11.9,12.65,15.48,17.11,18.05,18.58,19.05,19.12,19.01,18.81,18.56,18.28,17.99,16.52,15.18,14.01,12.99,11.32,1.02,8.974,8.125,7.421,6.829,6.325};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_12[i]*dens[j]*factor);}
G4double U_13[31]={5.596,6.24,7.447,8.566,9.609,1.58,11.5,12.35,13.15,16.27,18.17,19.29,19.97,2.61,2.78,2.74,2.59,2.37,2.11,19.83,18.36,16.97,15.75,14.67,12.87,11.46,1.31,9.371,8.585,7.92,7.35};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_13[i]*dens[j]*factor);}
G4double U_14[31]={5.788,6.446,7.681,8.832,9.908,1.92,11.87,12.76,13.6,16.98,19.13,2.44,21.26,22.08,22.37,22.4,22.29,22.1,21.87,21.6,2.14,18.72,17.45,16.32,14.42,12.9,11.66,1.63,9.77,9.036,8.404};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_14[i]*dens[j]*factor);}
G4double U_15[31]={6.009,6.682,7.945,9.126,1.23,11.28,12.26,13.19,14.07,17.7,2.09,21.6,22.57,23.6,24.02,24.14,24.09,23.95,23.75,23.51,22.09,2.66,19.34,18.16,16.14,14.51,13.17,12.05,11.1,1.29,9.584};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_15[i]*dens[j]*factor);}
G4double U_16[31]={6.191,6.883,8.18,9.394,1.54,11.61,12.63,13.6,14.52,18.37,21.01,22.72,23.85,25.1,25.65,25.87,25.89,25.8,25.64,25.42,24.07,22.63,21.28,2.06,17.94,16.2,14.76,13.54,12.51,11.61,1.84};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_16[i]*dens[j]*factor);}
G4double U_17[31]={6.421,7.126,8.445,9.68,1.85,11.95,13,13.99,14.94,18.98,21.85,23.76,25.04,26.5,27.19,27.5,27.59,27.55,27.42,27.24,25.94,24.51,23.14,21.87,19.66,17.83,16.3,15,13.89,12.92,12.08};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_17[i]*dens[j]*factor);}
G4double U_18[31]={6.564,7.283,8.624,9.878,11.06,12.19,13.26,14.28,15.25,19.47,22.54,24.64,26.08,27.75,28.59,28.99,29.15,29.15,29.06,28.9,27.66,26.23,24.83,23.53,21.25,19.35,17.74,16.37,15.19,14.17,13.27};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],U_18[i]*dens[j]*factor);}
j++;

G4double V_3[31]={2.541,2.772,3.148,3.426,3.623,3.759,3.849,3.907,3.942,3.94,3.827,3.688,3.546,3.276,3.037,2.827,2.642,2.479,2.335,2.206,1.728,1.422,1.21,1.055,.8447,.7079,.6118,.5404,.4851,.441,.4049};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_3[i]*dens[j]*factor);}
G4double V_4[31]={3.127,3.434,3.963,4.392,4.729,4.984,5.174,5.313,5.413,5.594,5.553,5.441,5.304,5.014,4.733,4.473,4.236,4.019,3.822,3.642,2.941,2.463,2.119,1.86,1.498,1.259,1.088,.9614,.8628,.784,.7195};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_4[i]*dens[j]*factor);}
G4double V_5[31]={3.591,3.966,4.629,5.192,5.661,6.041,6.342,6.577,6.759,7.191,7.263,7.208,7.1,6.825,6.532,6.246,5.974,5.72,5.483,5.263,4.366,3.719,3.236,2.863,2.328,1.964,1.702,1.504,1.35,1.227,1.126};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_5[i]*dens[j]*factor);}
G4double V_6[31]={3.96,4.395,5.175,5.856,6.445,6.946,7.362,7.702,7.977,8.714,8.935,8.961,8.899,8.665,8.38,8.085,7.795,7.517,7.251,7,5.94,5.14,4.523,4.035,3.315,2.814,2.446,2.166,1.946,1.769,1.623};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_6[i]*dens[j]*factor);}
G4double V_7[31]={4.304,4.787,5.659,6.433,7.122,7.727,8.249,8.691,9.061,1.14,1.55,1.68,1.68,1.51,1.26,9.969,9.674,9.383,9.099,8.826,7.635,6.699,5.955,5.353,4.446,3.799,3.316,2.944,2.649,2.41,2.213};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_7[i]*dens[j]*factor);}
G4double V_8[31]={4.602,5.127,6.083,6.94,7.716,8.414,9.032,9.572,1.04,11.48,12.09,12.35,12.44,12.36,12.15,11.89,11.6,11.31,11.01,1.73,9.436,8.381,7.52,6.81,5.715,4.916,4.311,3.839,3.461,3.152,2.896};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_8[i]*dens[j]*factor);}
G4double V_9[31]={4.84,5.401,6.429,7.356,8.203,8.977,9.679,1.3,1.85,12.66,13.5,13.88,14.05,14.08,13.92,13.68,13.41,13.12,12.82,12.53,11.16,1.02,9.063,8.264,7.01,6.077,5.359,4.792,4.334,3.958,3.643};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_9[i]*dens[j]*factor);}
G4double V_10[31]={5.057,5.649,6.739,7.728,8.637,9.476,1.25,1.95,11.58,13.75,14.82,15.35,15.62,15.76,15.66,15.46,15.2,14.92,14.62,14.33,12.91,11.68,1.65,9.768,8.365,7.304,6.477,5.817,5.278,4.832,4.457};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_10[i]*dens[j]*factor);}
G4double V_11[31]={5.234,5.85,6.992,8.036,8.999,9.896,1.73,11.5,12.21,14.75,16.08,16.8,17.18,17.46,17.44,17.29,17.06,16.8,16.51,16.22,14.77,13.48,12.36,11.4,9.853,8.658,7.716,6.955,6.33,5.809,5.367};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_11[i]*dens[j]*factor);}
G4double V_12[31]={5.449,6.079,7.253,8.337,9.343,1.28,11.16,11.99,12.75,15.62,17.22,18.12,18.62,19.06,19.13,19.03,18.84,18.59,18.32,18.04,16.57,15.22,14.05,13.02,11.34,1.03,8.983,8.13,7.423,6.83,6.324};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_12[i]*dens[j]*factor);}
G4double V_13[31]={5.647,6.293,7.503,8.623,9.666,1.64,11.57,12.43,13.25,16.42,18.29,19.37,2.01,2.63,2.8,2.76,2.61,2.4,2.15,19.87,18.41,17.02,15.79,14.7,12.9,11.47,1.32,9.377,8.588,7.921,7.35};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_13[i]*dens[j]*factor);}
G4double V_14[31]={5.843,6.505,7.745,8.897,9.973,1.98,11.94,12.84,13.7,17.14,19.26,2.53,21.31,22.1,22.38,22.42,22.31,22.13,21.9,21.64,2.19,18.77,17.49,16.35,14.44,12.91,11.67,1.64,9.772,9.036,8.403};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_14[i]*dens[j]*factor);}
G4double V_15[31]={6.072,6.75,8.021,9.205,1.31,11.35,12.34,13.28,14.17,17.87,2.25,21.71,22.64,23.63,24.04,24.16,24.12,23.98,23.79,23.55,22.15,2.71,19.39,18.2,16.17,14.53,13.18,12.06,11.1,1.29,9.583};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_15[i]*dens[j]*factor);}
G4double V_16[31]={6.259,6.957,8.265,9.483,1.63,11.7,12.72,13.69,14.62,18.55,21.18,22.85,23.93,25.13,25.68,25.89,25.92,25.84,25.68,25.47,24.14,22.7,21.35,2.11,17.98,16.23,14.78,13.55,12.51,11.62,1.84};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_16[i]*dens[j]*factor);}
G4double V_17[31]={6.496,7.208,8.54,9.782,1.95,12.05,13.09,14.09,15.04,19.17,22.04,23.91,25.14,26.54,27.22,27.53,27.62,27.59,27.47,27.29,26.01,24.58,23.2,21.92,19.7,17.86,16.32,15.01,13.89,12.93,12.08};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_17[i]*dens[j]*factor);}
G4double V_18[31]={6.642,7.369,8.724,9.986,11.17,12.3,13.36,14.38,15.35,19.64,22.74,24.8,26.18,27.8,28.61,29,29.17,29.18,29.09,28.94,27.73,26.29,24.89,23.58,21.29,19.37,17.75,16.37,15.19,14.16,13.27};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],V_18[i]*dens[j]*factor);}
j++;

G4double W_3[31]={1.603,1.75,1.999,2.196,2.352,2.474,2.568,2.64,2.694,2.799,2.774,2.705,2.623,2.454,2.295,2.151,2.023,1.908,1.804,1.711,1.358,1.127,.965,.8452,.6807,.573,.4969,.4401,.396,.3606,.3316};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_3[i]*dens[j]*factor);}
G4double W_4[31]={1.976,2.173,2.518,2.808,3.052,3.255,3.421,3.556,3.665,3.952,4.014,3.986,3.921,3.751,3.572,3.4,3.238,3.088,2.949,2.821,2.31,1.952,1.69,1.49,1.209,1.02,.8848,.7835,.7046,.6415,.5896};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_4[i]*dens[j]*factor);}
G4double W_5[31]={2.267,2.509,2.939,3.312,3.636,3.916,4.156,4.36,4.531,5.044,5.229,5.268,5.24,5.101,4.924,4.74,4.56,4.387,4.224,4.07,3.425,2.946,2.581,2.295,1.879,1.593,1.385,1.227,1.104,1.005,.923};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_5[i]*dens[j]*factor);}
G4double W_6[31]={2.498,2.779,3.287,3.735,4.133,4.486,4.798,5.069,5.305,6.07,6.406,6.535,6.56,6.474,6.317,6.136,5.949,5.764,5.584,5.412,4.66,4.073,3.609,3.237,2.679,2.285,1.994,1.77,1.594,1.451,1.333};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_6[i]*dens[j]*factor);}
G4double W_7[31]={2.696,3.008,3.579,4.091,4.554,4.973,5.349,5.685,5.983,7.012,7.524,7.764,7.86,7.853,7.732,7.566,7.384,7.196,7.008,6.824,5.991,5.311,4.757,4.3,3.6,3.092,2.709,2.411,2.174,1.981,1.821};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_7[i]*dens[j]*factor);}
G4double W_8[31]={2.869,3.207,3.834,4.403,4.924,5.401,5.838,6.234,6.591,7.882,8.584,8.951,9.133,9.225,9.155,9.017,8.848,8.665,8.477,8.289,7.399,6.641,6.006,5.472,4.632,4.007,3.527,3.149,2.845,2.595,2.387};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_8[i]*dens[j]*factor);}
G4double W_9[31]={3.02,3.38,4.053,4.67,5.241,5.769,6.256,6.705,7.116,8.654,9.546,1.05,1.32,1.52,1.5,1.4,1.25,1.07,9.886,9.696,8.768,7.95,7.249,6.65,5.69,4.961,4.393,3.939,3.57,3.265,3.009};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_9[i]*dens[j]*factor);}
G4double W_10[31]={3.16,3.539,4.25,4.909,5.523,6.095,6.629,7.125,7.584,9.359,1.44,11.09,11.46,11.78,11.83,11.76,11.63,11.47,11.29,11.1,1.15,9.283,8.527,7.87,6.8,5.972,5.317,4.789,4.355,3.993,3.688};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_10[i]*dens[j]*factor);}
G4double W_11[31]={3.269,3.662,4.406,5.102,5.754,6.367,6.942,7.482,7.987,9.996,11.28,12.08,12.57,13.04,13.17,13.15,13.05,12.91,12.74,12.56,11.6,1.7,9.889,9.18,8.003,7.076,6.334,5.727,5.225,4.803,4.444};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_11[i]*dens[j]*factor);}
G4double W_12[31]={3.392,3.795,4.561,5.284,5.967,6.612,7.222,7.798,8.341,1.56,12.04,12.99,13.6,14.22,14.44,14.47,14.41,14.29,14.14,13.96,13,12.07,11.22,1.47,9.203,8.191,7.369,6.692,6.126,5.647,5.237};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_12[i]*dens[j]*factor);}
G4double W_13[31]={3.506,3.919,4.705,5.45,6.159,6.833,7.473,8.081,8.659,11.06,12.73,13.85,14.58,15.37,15.69,15.79,15.77,15.68,15.54,15.38,14.44,13.49,12.61,11.81,1.46,9.366,8.467,7.719,7.088,6.551,6.089};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_13[i]*dens[j]*factor);}
G4double W_14[31]={3.618,4.04,4.843,5.607,6.337,7.034,7.701,8.337,8.945,11.52,13.37,14.63,15.49,16.46,16.89,17.05,17.08,17.02,16.91,16.76,15.84,14.88,13.97,13.14,11.71,1.54,9.573,8.76,8.07,7.478,6.966};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_14[i]*dens[j]*factor);}
G4double W_15[31]={3.732,4.163,4.982,5.763,6.513,7.233,7.924,8.587,9.221,11.96,13.98,15.41,16.4,17.56,18.11,18.37,18.45,18.43,18.36,18.24,17.38,16.41,15.48,14.63,13.12,11.87,1.82,9.936,9.178,8.524,7.954};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_15[i]*dens[j]*factor);}
G4double W_16[31]={3.835,4.276,5.116,5.916,6.686,7.428,8.142,8.829,9.488,12.38,14.56,16.15,17.28,18.64,19.33,19.67,19.82,19.85,19.81,19.72,18.92,17.97,17.03,16.15,14.58,13.26,12.13,11.17,1.35,9.631,9.005};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_16[i]*dens[j]*factor);}
G4double W_17[31]={3.958,4.407,5.258,6.072,6.857,7.616,8.349,9.056,9.738,12.76,15.1,16.83,18.09,19.65,2.48,2.91,21.12,21.2,21.2,21.13,2.41,19.47,18.52,17.62,15.99,14.6,13.41,12.39,11.5,1.73,1.05};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_17[i]*dens[j]*factor);}
G4double W_18[31]={4.048,4.504,5.369,6.194,6.991,7.764,8.512,9.236,9.936,13.07,15.55,17.43,18.81,2.57,21.52,22.05,22.33,22.45,22.48,22.44,21.78,2.85,19.89,18.97,17.29,15.85,14.6,13.53,12.59,11.77,11.05};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],W_18[i]*dens[j]*factor);}
j++;

G4double X_3[31]={3.332,3.612,4.056,4.373,4.592,4.739,4.835,4.893,4.924,4.88,4.713,4.519,4.326,3.97,3.659,3.39,3.157,2.952,2.773,2.613,2.031,1.663,1.411,1.228,.9798,.8194,.7071,.6238,.5594,.508,.466};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_3[i]*dens[j]*factor);}
G4double X_4[31]={4.135,4.515,5.152,5.649,6.03,6.313,6.519,6.667,6.771,6.928,6.834,6.663,6.468,6.07,5.698,5.36,5.055,4.781,4.533,4.309,3.45,2.875,2.465,2.159,1.735,1.455,1.256,1.109,.9941,.9026,.8279};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_4[i]*dens[j]*factor);}
G4double X_5[31]={4.766,5.239,6.053,6.721,7.261,7.69,8.024,8.281,8.476,8.912,8.94,8.827,8.657,8.263,7.862,7.483,7.129,6.802,6.501,6.223,5.115,4.335,3.759,3.318,2.69,2.266,1.961,1.732,1.554,1.411,1.294};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_5[i]*dens[j]*factor);}
G4double X_6[31]={5.255,5.811,6.786,7.61,8.302,8.877,9.345,9.723,1.02,1.8,1.99,1.97,1.84,1.48,1.08,9.68,9.295,8.932,8.59,8.271,6.952,5.981,5.244,4.666,3.822,3.239,2.813,2.489,2.236,2.031,1.863};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_6[i]*dens[j]*factor);}
G4double X_7[31]={5.72,6.34,7.441,8.391,9.212,9.916,1.51,11.01,11.42,12.58,12.97,13.06,13,12.71,12.33,11.93,11.53,11.14,1.77,1.42,8.924,7.783,6.891,6.178,5.115,4.363,3.805,3.376,3.037,2.763,2.536};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_7[i]*dens[j]*factor);}
G4double X_8[31]={6.114,6.792,8.008,9.073,1.01,1.83,11.55,12.16,12.68,14.26,14.89,15.11,15.14,14.95,14.61,14.22,13.82,13.42,13.03,12.66,11.02,9.728,8.693,7.849,6.564,5.636,4.937,4.394,3.96,3.606,3.313};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_8[i]*dens[j]*factor);}
G4double X_9[31]={6.408,7.134,8.448,9.611,1.65,11.58,12.4,13.11,13.74,15.73,16.61,16.98,17.1,17.01,16.72,16.36,15.96,15.56,15.17,14.78,13.04,11.63,1.47,9.522,8.047,6.961,6.132,5.48,4.955,4.523,4.163};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_9[i]*dens[j]*factor);}
G4double X_10[31]={6.674,7.439,8.836,1.09,11.21,12.23,13.14,13.96,14.68,17.09,18.23,18.77,19,19.03,18.81,18.47,18.09,17.7,17.3,16.9,15.08,13.56,12.31,11.26,9.601,8.365,7.408,6.648,6.03,5.519,5.09};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_10[i]*dens[j]*factor);}
G4double X_11[31]={6.899,7.695,9.16,1.48,11.69,12.78,13.78,14.69,15.51,18.36,19.81,2.54,2.91,21.1,2.96,2.68,2.33,19.94,19.55,19.15,17.27,15.66,14.3,13.15,11.31,9.914,8.822,7.945,7.227,6.629,6.124};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_11[i]*dens[j]*factor);}
G4double X_12[31]={7.191,8.003,9.507,1.88,12.14,13.3,14.36,15.34,16.23,19.48,21.23,22.17,22.68,23.05,23.01,22.78,22.46,22.1,21.71,21.32,19.39,17.71,16.26,15.03,13.03,11.49,1.27,9.286,8.473,7.791,7.213};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_12[i]*dens[j]*factor);}
G4double X_13[31]={7.454,8.287,9.835,11.26,12.57,13.78,14.9,15.94,16.9,2.51,22.57,23.72,24.38,24.95,25.02,24.86,24.59,24.26,23.89,23.5,21.56,19.81,18.29,16.97,14.82,13.14,11.8,1.71,9.799,9.032,8.377};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_13[i]*dens[j]*factor);}
G4double X_14[31]={7.717,8.569,1.15,11.61,12.97,14.22,15.39,16.49,17.5,21.44,23.79,25.16,25.97,26.73,26.92,26.84,26.62,26.32,25.97,25.6,23.65,21.85,2.27,18.88,16.6,14.8,13.34,12.15,11.15,1.3,9.576};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_14[i]*dens[j]*factor);}
G4double X_15[31]={8.04,8.913,1.54,12.04,13.43,14.73,15.94,17.08,18.14,22.39,25.04,26.63,27.61,28.59,28.92,28.93,28.78,28.52,28.21,27.86,25.94,24.11,22.46,21.01,18.58,16.64,15.06,13.75,12.65,11.71,1.91};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_15[i]*dens[j]*factor);}
G4double X_16[31]={8.294,9.193,1.86,12.4,13.84,15.18,16.44,17.63,18.74,23.29,26.23,28.06,29.21,3.42,3.9,31.02,3.94,3.74,3.46,3.13,28.27,26.41,24.72,23.2,2.63,18.56,16.86,15.44,14.24,13.21,12.32};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_16[i]*dens[j]*factor);}
G4double X_17[31]={8.625,9.54,11.24,12.8,14.27,15.64,16.94,18.16,19.31,24.1,27.32,29.37,3.69,32.13,32.76,32.97,32.96,32.81,32.57,32.27,3.46,28.59,26.85,25.28,22.6,2.42,18.61,17.09,15.8,14.69,13.72};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_17[i]*dens[j]*factor);}
G4double X_18[31]={8.813,9.748,11.47,13.07,14.56,15.96,17.28,18.53,19.71,24.73,28.21,3.48,31.98,33.65,34.43,34.74,34.8,34.7,34.49,34.23,32.47,3.59,28.82,27.21,24.43,22.15,2.26,18.65,17.29,16.1,15.07};
pv=new G4LPhysicsFreeVector(31,E[0],E[30]);pv->SetSpline(spline);dedx.push_back(pv);
for(i=0;i<31;i++){pv->PutValues(i,E[i],X_18[i]*dens[j]*factor);}
j++;



if(corr) {
  G4int n = 0;
  for(j=0; j<31; j++) {
    for(i=0; i<16; i++) {
      corr->AddStoppingData(Z[i], AA_Ion[i], NameMaterial[j], dedx[n]);
      n++;
    }
  }
}
}
