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
// $Id: G4WaterStoppingRange.cc,v 1.2 2008-08-14 16:05:18 antoni Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Description: Data on stopping power
//
// Author:      Ivanchenko 1.07.2008
//
// Modifications:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4WaterStoppingRange.hh"
#include "G4LPhysicsFreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4WaterStoppingRange::G4WaterStoppingRange(G4bool splineFlag) 
{
  spline = splineFlag;
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4WaterStoppingRange::~G4WaterStoppingRange()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4WaterStoppingRange::GetRange(G4int ionZ, G4double kinEnergy)
{
  G4double res = 0.0;
  if(ionZ < 3 || ionZ > 18 ) return res; 
  G4bool b;
  G4int idx = ionZ - 3;
  G4cout << "Index  " << idx << G4endl;
  G4cout << "Number " << R.size() << G4endl;
  G4cout << "Pointer  " << R[idx] << G4endl;
  G4double scaledEnergy = kinEnergy/A[ionZ - 3];
  G4cout << "Scaled Energy is " << scaledEnergy << G4endl;
  G4double emin = 0.025*MeV;
  if(scaledEnergy < emin) {
    res = (R[idx])->GetValue(emin, b)*sqrt(scaledEnergy/emin);
  } else {
    res = (R[idx])->GetValue(scaledEnergy, b);
  }
  return res;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PhysicsVector* 
G4WaterStoppingRange::GetPhysicsVector(G4int ionZ)
{
  if(ionZ < 3 || ionZ > 18) return 0; 
  G4int idx = ionZ - 3;
  return R[idx];
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WaterStoppingRange::Initialise()
{
G4int i;
R.reserve(16*1);

//..List of ions
G4double factor = cm/1000;
G4LPhysicsFreeVector* pv;  
G4int Z_Ion[16] = {3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
G4double A_Ion[16] = {6.941,9.0122,10.811,12.011,14.007,15.999,18.998,20.180,22.990,24.305,26.982,28.086,30.974,32.065,35.453,39.948};

for(i=0; i<16; i++) { 
  Z[i] = Z_Ion[i];
  A[i] = A_Ion[i];
}

  //..Reduced energies
G4double E[53] = {.025,.03,.04,.05,.06,.07,.08,.09,.1,.15,.2,.25,.3,.4,.5,.6,.7,.8,.9,1,1.5,2,2.5,3,4,5,6,7,8,9,10,15,20,25,30,40,50,60,70,80,90,100,150,200,250,300,400,500,600,700,800,900,1000};
for(i=0; i<53; i++) {E[i] *= MeV;}

G4double Water_Li[53]={0,.01283,.03612,.0572,.07686,.09558,.1137,.1313,.1487,.2345,.3209,.4101,.5028,.6999,.9133,1.143,1.39,1.655,1.936,2.234,3.984,6.16,8.758,11.77,19.01,27.83,38.17,49.99,63.26,77.95,94.02,194.2,326.1,488.1,678.9,1142,1710,2376,3136,3985,4919,5935,12110,19930,29130,39520,63130,89820,118900,149700,182000,215500,249800};
pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
pv->SetSpline(spline);R.push_back(pv);
for(i=0; i<53; i++) {pv->PutValues(i,E[i],Water_Li[i]*factor);}

G4double Water_Be[53]={0,.01319,.03689,.05806,.07753,.09583,.1133,.1302,.1466,.226,.3036,.3821,.4625,.6302,.8081,.9968,1.197,1.408,1.63,1.864,3.207,4.839,6.759,8.968,14.24,2.63,28.11,36.66,46.26,56.88,68.51,141,236.5,353.8,491.9,827.3,1238,1720,2269,2883,3559,4294,8759,14410,21060,28570,45630,64920,85910,108200,131600,155700,180500};
pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
pv->SetSpline(spline);R.push_back(pv);
for(i=0;i<53;i++) {pv->PutValues(i,E[i],Water_Be[i]*factor);}

G4double Water_B[53]={0,.01393,.03872,.06061,.08054,.0991,.1166,.1334,.1496,.2265,.3,.373,.4469,.5988,.7576,.9239,1.098,1.28,1.471,1.67,2.791,4.125,5.673,7.436,11.61,16.63,22.49,29.18,36.68,44.99,54.07,11.8,185.5,277.3,385.5,648,969.3,1346,1776,2256,2784,3359,6851,11270,16470,22340,35670,50750,67160,84590,102800,121700,141100};
pv =new G4LPhysicsFreeVector(53,E[0],E[52]);
pv->SetSpline(spline);R.push_back(pv);
for(i=0; i<53; i++){pv->PutValues(i,E[i],Water_B[i]*factor);}

G4double Water_C[53]={0,.01375,.03802,.0592,.07832,.09597,.1125,.1282,.1433,.2138,.2797,.3443,.4089,.5405,.6763,.8172,.9636,1.116,1.273,1.437,2.346,3.408,4.624,5.997,9.21,13.05,17.51,22.58,28.27,34.55,41.43,84.36,141,21.6,292.6,491.6,735.1,1021,1346,1710,2110,2546,5189,8535,12470,16910,27000,38410,50830,64020,77830,92110,106800};
pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
pv->SetSpline(spline);R.push_back(pv);
for(i=0; i<53; i++) {pv->PutValues(i,E[i],Water_C[i]*factor);}

G4double Water_N[53]={0,.01482,.04081,.06333,.08353,.1021,.1194,.1357,.1513,.2231,.2889,.3527,.4159,.5433,.6735,.8074,.9455,1.088,1.235,1.387,2.221,3.179,4.263,5.475,8.287,11.62,15.47,19.83,24.71,3.1,35.99,72.73,121.2,18.9,251.1,421.8,63.7,875.7,1155,1467,1810,2183,4449,7315,10690,14490,23140,32910,43540,54840,66660,78890,91450};
pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
pv->SetSpline(spline);R.push_back(pv);
for(i=0; i<53; i++) {pv->PutValues(i,E[i],Water_N[i]*factor);}

G4double Water_O[53]={0,.01589,.04363,.06751,.08882,.1083,.1263,.1433,.1594,.2329,.2991,.3625,.4248,.5493,.6754,.8043,.9364,1.072,1.212,1.355,2.132,3.013,4,5.094,7.609,1.56,13.96,17.79,22.07,26.78,31.93,64.01,106.4,158.5,220,369.3,552.1,766.6,1011,1284,1584,1911,3893,6400,9349,12680,20240,28780,38080,47960,58290,68980,79950};
pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
pv->SetSpline(spline);R.push_back(pv);
for(i=0; i<53; i++) {pv->PutValues(i,E[i],Water_O[i]*factor);}

G4double Water_F[53]={0,.01804,.04939,.07624,.1001,.1217,.1417,.1605,.1782,.2582,.3292,.3965,.4623,.5926,.7237,.8569,.9928,1.132,1.274,1.42,2.204,3.083,4.057,5.131,7.575,1.42,13.67,17.32,21.38,25.84,3.7,6.9,10.7,149.6,207.3,347.5,519.1,72.3,949.7,1206,1488,1794,3653,6005,8771,11890,18980,26990,35710,44970,54660,64680,74970};
pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
pv->SetSpline(spline);R.push_back(pv);
for(i=0; i<53; i++) {pv->PutValues(i,E[i],Water_F[i]*factor);}

G4double Water_Ne[53]={0,.01827,.04995,.07696,.1008,.1224,.1423,.1609,.1784,.2568,.3254,.3898,.4523,.5754,.6983,.8226,.9489,1.078,1.209,1.343,2.06,2.855,3.731,4.689,6.853,9.35,12.18,15.36,18.87,22.71,26.9,52.79,86.77,128.5,177.7,297.3,443.7,615.4,811,1030,1270,1531,3116,5121,7478,10140,16180,23000,30430,38320,46570,55100,63870};
pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
pv->SetSpline(spline);R.push_back(pv);
for(i=0; i<53; i++) {pv->PutValues(i,E[i],Water_Ne[i]*factor);}

G4double Water_Na[53]={0,.02033,.05547,.08529,.1115,.1351,.1568,.177,.1959,.2799,.3524,.4196,.4842,.6104,.7351,.8604,.987,1.115,1.246,1.379,2.08,2.849,3.687,4.598,6.641,8.983,11.63,14.58,17.84,21.42,25.31,49.39,81.14,12.3,166.6,279.6,418.4,581.3,767.1,974.7,1203,1451,2958,4862,7101,9626,15360,21840,28890,36380,44220,52330,60650};
pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
pv->SetSpline(spline);R.push_back(pv);
for(i=0; i<53; i++) { pv->PutValues(i,E[i],Water_Na[i]*factor); }

G4double Water_Mg[53]={0,.02051,.05606,.08627,.1128,.1367,.1586,.1789,.198,.2821,.354,.4202,.4835,.6065,.7278,.8492,.9717,1.096,1.222,1.35,2.025,2.762,3.563,4.428,6.357,8.553,11.02,13.75,16.76,2.04,23.58,45.35,73.73,108.5,149.5,249.1,371.2,514.3,677.4,859.6,1060,1278,2598,4268,6230,8442,13470,19150,25320,31880,38750,45850,53140};
pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
pv->SetSpline(spline);R.push_back(pv);
for(i=0; i<53; i++) {pv->PutValues(i,E[i],Water_Mg[i]*factor);}

G4double Water_Al[53]={0,.02236,.06113,.09409,.123,.149,.1728,.1948,.2154,.3059,.3826,.4525,.519,.6475,.7732,.8986,1.025,1.152,1.281,1.412,2.099,2.842,3.645,4.51,6.425,8.59,11.01,13.68,16.61,19.79,23.23,44.24,71.52,104.9,144.2,239.7,356.6,493.8,65.2,824.9,1017,1226,2492,4091,5972,8091,12910,18340,24260,30540,37120,43920,50900};
pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
pv->SetSpline(spline);R.push_back(pv);
for(i=0; i<53; i++) {pv->PutValues(i,E[i],Water_Al[i]*factor);}

G4double Water_Si[53]={0,.02249,.06155,.09477,.1239,.1501,.174,.1961,.2168,.307,.3829,.4516,.5166,.6413,.7628,.8834,1.004,1.126,1.249,1.374,2.025,2.725,3.478,4.285,6.064,8.063,1.29,12.73,15.41,18.3,21.42,4.4,64.92,94.85,130,215.5,32.1,442.8,582.7,738.9,91.6,1097,2229,3659,5339,7234,11540,16400,21680,27300,33170,39240,45480};
pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
pv->SetSpline(spline);R.push_back(pv);
for(i=0; i<53; i++) {pv->PutValues(i,E[i],Water_Si[i]*factor);}

G4double Water_P[53]={0,.02404,.0659,.1016,.1329,.161,.1866,.2103,.2324,.3287,.409,.4813,.5493,.6791,.8047,.9288,1.053,1.178,1.303,1.43,2.09,2.795,3.549,4.354,6.118,8.092,1.28,12.67,15.28,18.1,21.13,39.51,63.17,92.02,125.9,208.2,309,427.3,562.2,712.7,878.4,1058,2150,3529,5149,6975,11120,15800,20900,26310,31970,37820,43830};
pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
pv->SetSpline(spline);R.push_back(pv);
for(i=0; i<53; i++) {pv->PutValues(i,E[i],Water_P[i]*factor);}

G4double Water_S[53]={0,.02414,.06619,.102,.1335,.1617,.1874,.2112,.2333,.3294,.409,.4803,.5469,.6733,.795,.9148,1.034,1.154,1.274,1.395,2.023,2.689,3.398,4.152,5.797,7.626,9.643,11.85,14.24,16.82,19.6,36.33,57.8,83.94,114.6,189.2,28.5,387.7,509.9,646.5,796.7,960,1950,3201,4670,6326,10090,14330,18950,23850,28980,34290,39730};
pv = new G4LPhysicsFreeVector(53,E[0],E[52]);pv->SetSpline(spline);R.push_back(pv);
for(i=0; i<53; i++) {pv->PutValues(i,E[i],Water_S[i]*factor);}

G4double Water_Cl[53]={0,.02551,.07006,.1082,.1416,.1717,.1991,.2243,.2479,.3498,.4338,.5087,.5783,.7097,.8355,.959,1.082,1.204,1.327,1.451,2.09,2.765,3.48,4.237,5.882,7.704,9.704,11.88,14.24,16.79,19.51,35.86,56.75,82.13,111.9,184.1,272.6,376.4,494.8,627.1,772.6,93.7,1890,3101,4524,6127,9767,13880,18350,23090,28060,33190,38460};
pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
pv->SetSpline(spline);R.push_back(pv);
for(i=0; i<53; i++) {pv->PutValues(i,E[i],Water_Cl[i]*factor);}

G4double Water_Ar[53]={0,.02863,.07865,.1215,.1591,.1929,.2236,.252,.2784,.3924,.4859,.5688,.6454,.7895,.9267,1.061,1.194,1.327,1.46,1.593,2.28,3.003,3.766,4.573,6.318,8.242,1.35,12.63,15.1,17.76,2.6,37.55,59.08,85.16,115.7,189.5,279.9,385.8,506.6,641.5,789.9,951.2,1929,3164,4614,6248,9956,14140,18690,23530,28580,33810,39180};
pv = new G4LPhysicsFreeVector(53,E[0],E[52]);
pv->SetSpline(spline);R.push_back(pv);
for(i=0; i<53; i++) {pv->PutValues(i,E[i],Water_Ar[i]*factor);}

}
