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
// Geant4 class G4NuclearRadii
//
// Author V.Ivanchenko 27.05.2019
//
//

#include "G4NuclearRadii.hh"
#include "G4Pow.hh"
#include "G4PhysicalConstants.hh"

G4double G4NuclearRadii::ExplicitRadius(G4int Z, G4int A)
{
  G4double R = 0.0;
  // Special rms radii for light nucleii
  if(Z <= 4) {
    if(A == 1)                { R = 0.89*CLHEP::fermi; }// p
    else if(A == 2)           { R = 2.13*CLHEP::fermi; }// d
    else if(Z == 1 && A == 3) { R = 1.80*CLHEP::fermi; }// t
    else if(Z == 2 && A == 3) { R = 1.96*CLHEP::fermi; }// He3
    else if(Z == 2 && A == 4) { R = 1.68*CLHEP::fermi; }// He4
    else if(Z == 3)           { R = 2.40*CLHEP::fermi; }// Li7
    else if(Z == 4)           { R = 2.51*CLHEP::fermi; }// Be9
  }
  return R;
}

G4double G4NuclearRadii::Radius(G4int Z, G4int A)
{
  G4double R = ExplicitRadius(Z, A);
  if(0.0 == R) {
    if (A <= 50) {
      G4double y = 1.1;
      if( A <= 15)      { y = 1.26; }
      else if( A <= 20) { y = 1.19; }
      else if( A <= 30) { y = 1.12; }
      G4double x = G4Pow::GetInstance()->Z13(A);
      R = y*(x - 1./x);
    } else {
      R = G4Pow::GetInstance()->powZ(A, 0.27);
    }
    R *= CLHEP::fermi;
  }
  return R;
}

G4double G4NuclearRadii::RadiusRMS(G4int Z, G4int A)
{
  G4double R = ExplicitRadius(Z, A);
  if(0.0 == R) {
    R = 1.24*G4Pow::GetInstance()->powZ(A, 0.28)*CLHEP::fermi;
  }
  return R;
}

G4double G4NuclearRadii::RadiusNNGG(G4int Z, G4int A)
{
  G4double R = ExplicitRadius(Z, A);
  if(0.0 == R) {
    if(A > 20) {
      R = 1.08*G4Pow::GetInstance()->Z13(A)
	*(0.85 + 0.15*G4Exp(-(G4double)(A - 21)/40.));
    } else {
      R = 1.08*G4Pow::GetInstance()->Z13(A)
	*(1.0 + 0.3*G4Exp(-(G4double)(A - 21)/10.));
    }
    R *= CLHEP::fermi;
  }
  return R;
}

G4double G4NuclearRadii::RadiusHNGG(G4int A)
{
  G4double R = CLHEP::fermi;
  if(A > 20) {
    R *= 1.08*G4Pow::GetInstance()->Z13(A)
      *(0.8 + 0.2*G4Exp(-(G4double)(A - 20)/20.));
  } else {
    R *= 1.08*G4Pow::GetInstance()->Z13(A)
      *(1.0 + 0.1*G4Exp(-(G4double)(A - 20)/20.));
  }
  return R;
}

G4double G4NuclearRadii::RadiusKNGG(G4int A)
{
  return 1.3*CLHEP::fermi*G4Pow::GetInstance()->Z13(A);
}

G4double G4NuclearRadii::RadiusND(G4int A)
{
  G4double R = CLHEP::fermi;
  if(1 == A) { R *= 0.89; }
  else {
    G4double x = G4Pow::GetInstance()->Z13(A);
    if(A <= 3.) { x *= 0.8; }
    else { x *= 1.7; }
    R *= x;
  } 
  return R;
}

G4double G4NuclearRadii::RadiusCB(G4int Z, G4int A)
{
  G4double R = ExplicitRadius(Z, A);
  if(0.0 == R) {
    G4int z = std::min(Z, 92);
    R = r0[z]*G4Pow::GetInstance()->Z13(A)*CLHEP::fermi;
  }
  return R;
}

const G4double G4NuclearRadii::r0[] = {
 1.2,
 1.3, 1.3, 1.3, 1.3,1.17,1.54,1.65,1.71, 1.7,1.75, // 1-10
 1.7,1.57,1.53, 1.4, 1.3,1.30,1.44, 1.4, 1.4, 1.4, //11-20
 1.4, 1.4,1.46, 1.4, 1.4,1.46,1.55, 1.5,1.38,1.48, //21-30
 1.4, 1.4, 1.4,1.46, 1.4, 1.4, 1.4, 1.4, 1.4,1.45, //31-40
 1.4, 1.4, 1.4, 1.4, 1.4, 1.4,1.45,1.48, 1.4,1.52, //41-50
1.46, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.5, //51-60
 1.4, 1.4, 1.4, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.4, //61-70
 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3,1.33,1.43, //71-80
 1.3,1.32,1.34, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, //81-90
 1.3, 1.3};
