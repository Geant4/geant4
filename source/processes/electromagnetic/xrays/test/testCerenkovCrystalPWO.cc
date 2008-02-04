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
// 
//
//  
//
//  Test routine for cerenkov proximity geometry taking into account reflection
//
// History:
//
// 07.11.07, V. Grichine 

#include "G4ios.hh"
#include <fstream>
#include <cmath>
#include "globals.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"


#include <iomanip>

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4RegionStore.hh"
#include "G4MaterialTable.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"

#include "G4VXTRenergyLoss.hh"
#include "G4RegularXTRadiator.hh"
#include "G4TransparentRegXTRadiator.hh"
#include "G4GammaXTRadiator.hh"
#include "G4StrawTubeXTRadiator.hh"

#include "G4XTRGammaRadModel.hh"
#include "G4XTRRegularRadModel.hh"
#include "G4XTRTransparentRegRadModel.hh"

#include "G4SynchrotronRadiation.hh"
#include "G4Cerenkov.hh"

#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"

//////////////////////////////////////////////////////////////////
//
// Return refractive index of C6F14 vs. photon energy in eV

G4double GetRefractiveIndexA(G4double photonEnergy)
{
  G4double n;
  G4double a = 1.177;
  G4double b = 0.0172;
  n          = a + b*photonEnergy;
  return n;
}

//////////////////////////////////////////////////////////////////
//
// Return refractive index of quartz vs. photon energy in eV

G4double GetRefractiveIndexB(G4double photonEnergy)
{
  G4double result, n2;

  G4double E1 = 10.666;
  G4double E2 = 18.125;
  G4double F1 = 46.411;
  G4double F2 = 228.71;
  G4double E  = photonEnergy;

  n2 = 1. + F1/( E1*E1 - E*E ) + F2/( E2*E2 - E*E );

  result = std::sqrt(n2); 

  return result;
}

///////////////////////////////////////////////////////////////////
//
// Photon emitted in B transported to C

G4double TransmissionBC(G4double photonEnergy, G4double gamma)
{
  G4double beta, cosB, sinB, nB, thickB, cosC, sinC, tmp, rBC, result;

  if( gamma < 1.) gamma = 1.;
   
  beta = std::sqrt(1 - 1/gamma/gamma);
  nB   = GetRefractiveIndexB(photonEnergy);

  cosB = 1./nB/beta;
  
  if( cosB > 1.) cosB = 1.;

  sinB   = std::sqrt(1 - cosB*cosB);

  thickB = 0.5;  // in cm

  sinC   = nB*sinB;

  if( sinC > 1.) sinC = 1.;

  cosC   = std::sqrt(1 - sinC*sinC);

  tmp = ( nB*cosB - cosC )/( nB*cosB + cosC );

  rBC = tmp*tmp;

  if(rBC > 1.)  rBC = 1.;


  result = 369.77*thickB*2.;  // 2 eV of band

  result *= sinB*sinB;

  result *= 1. - rBC;

  return result;
}

///////////////////////////////////////////////////////////////////
//
// Photon emitted in A transported to C

G4double TransmissionABC(G4double photonEnergy, G4double gamma)
{
  G4double beta, cosA, sinA, cosB, sinB, nA, nB, thickA, cosC, sinC, tmp, rAB, rBC, result;

  if( gamma < 1.) gamma = 1.;
   
  beta = std::sqrt(1 - 1/gamma/gamma);
  nA   = GetRefractiveIndexA(photonEnergy);
  nB   = GetRefractiveIndexB(photonEnergy);

  cosA = 1./nA/beta;
  
  if( cosA > 1.) cosA = 1.;

  sinA   = std::sqrt(1 - cosA*cosA);

  thickA = 1.5;  // in cm

  sinB = nA*sinA/nB;  

  if( sinB > 1.) sinB = 1.;

  cosB   = std::sqrt(1 - sinB*sinB);

  tmp = ( nA*cosA - nB*cosB )/( nA*cosA + nB*cosB );

  rAB = tmp*tmp;

  if(rAB > 1.)  rAB = 1.;

  sinC   = nB*sinB;

  if( sinC > 1.) sinC = 1.;

  cosC   = std::sqrt(1 - sinC*sinC);

  tmp = ( nB*cosB - cosC )/( nB*cosB + cosC );

  rBC = tmp*tmp;

  if(rBC > 1.)  rBC = 1.;


  result = 369.77*thickA*2.;  // 2 eV of band

  result *= sinB*sinB;

  result *= 1. - rAB;

  result *= 1. - rBC;

  return result;
}


int main()
{
  G4int i, iMax = 20;
  G4double energy, refA, refB, gamma, numberBC, numberABC, beta, protonMass, protonMom;

  /*
  for(i=0;i<iMax;i++)
  {
    energy = 1. + i*8./iMax;
    refA = GetRefractiveIndexA(energy);
    refB = GetRefractiveIndexB(energy);
    G4cout<<"photon energy = "<<energy<<" eV; C6F14 = "<<refA<<"; quartz = "<<refB<<G4endl;  
  }
  

  energy = 7.;
  protonMass = 0.938; 


  for(i=0;i<iMax;i++)
  {
    gamma = 1.1 + i*2./iMax;
    beta = std::sqrt(1 - 1/gamma/gamma);
    protonMom = protonMass*beta*gamma;
    numberBC = TransmissionBC(energy, gamma);
    // G4cout<<"gamma = "<<gamma<<"; numberBC = "<<number<<G4endl;  
    G4cout<<"protonMom = "<<protonMom<<"; numberBC = "<<numberBC<<G4endl;  
  }

  for(i=0;i<iMax;i++)
  {
    gamma = 1.1 + i*2./iMax;
    beta = std::sqrt(1 - 1/gamma/gamma);
    protonMom = protonMass*beta*gamma;
    numberABC = TransmissionABC(energy, gamma);
    numberBC = TransmissionBC(energy, gamma);
    // G4cout<<"gamma = "<<gamma<<"; number = "<<number<<G4endl;  
    G4cout<<"protonMom = "<<protonMom<<"; numberBC = "<<numberBC
          <<"; numberABC = "<<numberABC<<G4endl;  
  }
  */

  // Arrays of CMS crystal dimensions, table 3.2 from ECAL TDR

  G4double cmsAF[17] = { 21.83, 21.83, 21.83, 21.83, 21.83, 
                         21.83, 21.83, 21.83, 21.83, 
                         21.83, 21.83, 21.83, 21.83,
                         21.83, 21.83, 21.83, 21.83         };

  G4double cmsBF[17] = { 23.59, 22.22, 22.34, 22.47, 22.61, 
                         22.60, 22.55, 22.67, 22.82,
                         23.08, 23.14, 23.29, 23.47, 
                         23.71, 23.88, 24.06, 24.29         };

  G4double cmsCF[17] = { 21.85, 21.87, 21.91, 21.94, 21.97, 
                         22.00, 22.03, 22.05, 22.08, 
                         22.10, 22.12, 22.14, 22.15, 
                         22.17, 22.18, 22.20, 22.21         };


  G4double cmsAR[17] = { 25.84, 25.81, 25.75, 25.67, 25.56,
                         25.43, 25.29, 25.14, 24.98, 
                         24.82, 24.65, 24.49, 24.33, 
                         24.17, 24.02, 23.88, 23.74         };

  G4double cmsBR[17] = { 25.48, 26.22, 26.28, 26.32, 26.34, 
                         26.18, 25.96, 25.92, 25.90, 
                         26.00, 25.89, 25.86, 25.87, 
                         25.95, 25.96, 25.99, 26.07         };

  G4double cmsCR[17] = { 25.86, 25.86, 25.84, 25.80, 25.72, 
                         25.63, 25.52, 25.39, 25.26, 
                         25.12, 24.97, 24.83, 24.68, 
                         24.54, 24.40, 24.27, 24.15         };

  iMax = 17;

  G4double meanAF = 0.,  meanBF = 0.,  meanCF = 0.,  
           meanAR = 0.,  meanBR = 0.,  meanCR = 0.;
  
  for(i=0;i<iMax;i++)
  {
    meanAF += cmsAF[i];
    meanBF += cmsBF[i];
    meanCF += cmsCF[i];
    meanAR += cmsAR[i];
    meanBR += cmsBR[i];
    meanCR += cmsCR[i];
  }
  meanAF /= iMax;
  meanCF /= iMax;
  meanBF /= iMax;
  meanAR /= iMax;
  meanBR /= iMax;
  meanCR /= iMax;

  G4double meanBFR  = 0.5*(meanBF + meanBR);
  G4double meanAFR  = 0.5*(meanAF + meanAR);
  G4double meanCFR  = 0.5*(meanCF + meanCR);
  G4double meanACFR = 0.5*(meanCFR + meanAFR);

  G4cout <<"meanAF = " << meanAF << "; meanBF = "<<meanBF << "; meanCF = " << meanCF << G4endl;

  G4cout << "meanAR = " << meanAR << "; meanBR = "<<meanBR << "; meanCR = "<<meanCR << G4endl;

  G4cout << "meanACFR = " << meanACFR << "; meanBFR = "<<meanBFR << G4endl;

  return 1 ;
}

