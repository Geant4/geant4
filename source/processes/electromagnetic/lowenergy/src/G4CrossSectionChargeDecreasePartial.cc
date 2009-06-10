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
// $Id: G4CrossSectionChargeDecreasePartial.cc,v 1.3 2009-06-10 13:32:36 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4CrossSectionChargeDecreasePartial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionChargeDecreasePartial::G4CrossSectionChargeDecreasePartial()
{

  //PROTON
  f0[0][0]=1.;
  a0[0][0]=-0.180;
  a1[0][0]=-3.600;
  b0[0][0]=-18.22;
  b1[0][0]=-1.997;
  c0[0][0]=0.215;
  d0[0][0]=3.550;
  x0[0][0]=3.450;
  x1[0][0]=5.251;

  numberOfPartialCrossSections[0] = 1;

  //ALPHA++
  f0[0][1]=1.;  a0[0][1]=0.95;
  a1[0][1]=-2.75;
  b0[0][1]=-23.00;
  c0[0][1]=0.215;
  d0[0][1]=2.95;
  x0[0][1]=3.50;

  f0[1][1]=1.;
  a0[1][1]=0.95;
  a1[1][1]=-2.75;
  b0[1][1]=-23.73;
  c0[1][1]=0.250;
  d0[1][1]=3.55;
  x0[1][1]=3.72;

  x1[0][1]=-1.;
  b1[0][1]=-1.;

  x1[1][1]=-1.;
  b1[1][1]=-1.;  
 
  numberOfPartialCrossSections[1] = 2;

  // ALPHA+
  f0[0][2]=1.;
  a0[0][2]=0.65;
  a1[0][2]=-2.75;
  b0[0][2]=-21.81;
  c0[0][2]=0.232;
  d0[0][2]=2.95;
  x0[0][2]=3.53;

  x1[0][2]=-1.;
  b1[0][2]=-1.;
 
  numberOfPartialCrossSections[2] = 1;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrossSectionChargeDecreasePartial::~G4CrossSectionChargeDecreasePartial()
{}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CrossSectionChargeDecreasePartial::CrossSection(G4double k, G4int index, 
							   const G4ParticleDefinition* particleDefinition)
{
  G4int particleTypeIndex = 0;
  G4DNAGenericIonsManager* instance;
  instance = G4DNAGenericIonsManager::Instance();

  if (particleDefinition == G4Proton::ProtonDefinition()) particleTypeIndex=0;
  if (particleDefinition == instance->GetIon("alpha++")) particleTypeIndex=1;
  if (particleDefinition == instance->GetIon("alpha+")) particleTypeIndex=2;

  //
  // sigma(T) = f0 10 ^ y(log10(T/eV))
  //
  //         /  a0 x + b0                    x < x0
  //         |
  // y(x) = <   a0 x + b0 - c0 (x - x0)^d0   x0 <= x < x1
  //         |
  //         \  a1 x + b1                    x >= x1
  //
  //
  // f0, a0, a1, b0, b1, c0, d0, x0, x1 are parameters that change for protons and helium (0, +, ++)
  //
  // f0 has been added to the code in order to manage partial (shell-dependent) cross sections (if no shell dependence is present. f0=1. Sum of f0 over the considered shells should give 1)
  //
  // From Rad. Phys. and Chem. 59 (2000) 255-275, M. Dingfelder et al.
  // Inelastic-collision cross sections of liquid water for interactions of energetic proton
  //
  
  if (x1[index][particleTypeIndex]<x0[index][particleTypeIndex])
  {
      //
      // if x1 < x0 means that x1 and b1 will be calculated with the following formula (this piece of code is run on all alphas and not on protons)
      //
      // x1 = x0 + ((a0 - a1)/(c0 * d0)) ^ (1 / (d0 - 1))
      //
      // b1 = (a0 - a1) * x1 + b0 - c0 * (x1 - x0) ^ d0
      //
 
      x1[index][particleTypeIndex]=x0[index][particleTypeIndex] + std::pow((a0[index][particleTypeIndex] - a1[index][particleTypeIndex]) 
									   / (c0[index][particleTypeIndex] * d0[index][particleTypeIndex]), 1. / (d0[index][particleTypeIndex] - 1.));
      b1[index][particleTypeIndex]=(a0[index][particleTypeIndex] - a1[index][particleTypeIndex]) * x1[index][particleTypeIndex] 
	+ b0[index][particleTypeIndex] - c0[index][particleTypeIndex] * std::pow(x1[index][particleTypeIndex] 
										 - x0[index][particleTypeIndex], d0[index][particleTypeIndex]);
  }

  G4double x(std::log10(k/eV));
  G4double y;
  
  if (x<x0[index][particleTypeIndex])
    y=a0[index][particleTypeIndex] * x + b0[index][particleTypeIndex];
  else if (x<x1[index][particleTypeIndex])
    y=a0[index][particleTypeIndex] * x + b0[index][particleTypeIndex] - c0[index][particleTypeIndex] 
      * std::pow(x - x0[index][particleTypeIndex], d0[index][particleTypeIndex]);
  else
    y=a1[index][particleTypeIndex] * x + b1[index][particleTypeIndex];

  return f0[index][particleTypeIndex] * std::pow(10., y)*m*m;
  
}

G4int G4CrossSectionChargeDecreasePartial::RandomSelect(G4double k, 
							const G4ParticleDefinition* particleDefinition)
{
  G4int particleTypeIndex = 0;
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  if (particleDefinition == G4Proton::ProtonDefinition()) particleTypeIndex = 0;
  if (particleDefinition == instance->GetIon("alpha++")) particleTypeIndex = 1;
  if (particleDefinition == instance->GetIon("alpha+")) particleTypeIndex = 2;

  const G4int n = numberOfPartialCrossSections[particleTypeIndex];
  G4double* values(new G4double[n]);
  G4double value(0);
  G4int i = n;
  
  while (i>0)
  {
      i--;
      values[i]=CrossSection(k, i, particleDefinition);
      value+=values[i];
  }
  
  value*=G4UniformRand();
  
  i=n;
  while (i>0)
  {
      i--;
   
      if (values[i]>value)
	break;
  
      value-=values[i];
  }
  
  delete[] values;
  
  return i;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CrossSectionChargeDecreasePartial::Sum(G4double k, const G4ParticleDefinition* particleDefinition)
{
  G4int particleTypeIndex = 0;
  G4DNAGenericIonsManager* instance;
  instance = G4DNAGenericIonsManager::Instance();

  if (particleDefinition == G4Proton::ProtonDefinition()) particleTypeIndex=0;
  if (particleDefinition == instance->GetIon("alpha++")) particleTypeIndex=1;
  if (particleDefinition == instance->GetIon("alpha+")) particleTypeIndex=2;

  G4double totalCrossSection = 0.;

  for (G4int i=0; i<numberOfPartialCrossSections[particleTypeIndex]; i++)
  {
    totalCrossSection += CrossSection(k,i,particleDefinition);
  }
  return totalCrossSection;
}
