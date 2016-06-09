//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4Generator2BS
//
// Author:        Andreia Trindade (andreia@lip.pt)
//                Pedro Rodrigues  (psilva@lip.pt)
//                Luis Peralta     (luis@lip.pt)
//
// Creation date: 2 June 2003
//
// Modifications: 
// 02 Jun 2003                    First implementation acording with new design
// 05 Nov 2003  MGP               Fixed std namespace
// 17 Nov 2003  MGP               Fixed compilation problem on Windows                  
//
// Class Description: 
//
// Concrete base class for Bremsstrahlung Angular Distribution Generation - 2BS Distribution
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//    

#include "G4Generator2BS.hh"
#include "Randomize.hh"
//    

G4Generator2BS::G4Generator2BS(const G4String& name):G4VBremAngularDistribution(name)
{;}

//    

G4Generator2BS::~G4Generator2BS() 
{;}

//

G4double G4Generator2BS::PolarAngle(const G4double initial_energy,
				    const G4double final_energy,
				    const G4int Z)
{

  // Adapted from "Improved bremsstrahlung photon angular sampling in the EGS4 code system"
  // by Alex F. Bielajew, Rahde Mohan anc Chen-Shou Chui, PIRS-0203
  // Ionizing Radiation Standards
  // Institute for National Measurement Standards 
  // National Research Council of Canada
  // Departement of Medical Physics, Memorial Sloan-Kettering Cancer Center, New York


  G4double theta = 0;

  G4double initialTotalEnergy = (initial_energy+electron_mass_c2)/electron_mass_c2;
  G4double finalTotalEnergy = (final_energy+electron_mass_c2)/electron_mass_c2;
  EnergyRatio = finalTotalEnergy/initialTotalEnergy;
  G4double gMaxEnergy = (M_PI*initialTotalEnergy)*(M_PI*initialTotalEnergy);

  G4double Zeff = sqrt(static_cast<G4double>(Z) * (static_cast<G4double>(Z) + 1.0));
  z = (0.00008116224*(pow(Zeff,0.3333333)));

  // Rejection arguments
  rejection_argument1 = (1.0+EnergyRatio*EnergyRatio); 
  rejection_argument2 = - 2.0*EnergyRatio + 3.0*rejection_argument1;
  rejection_argument3 = ((1-EnergyRatio)/(2.0*initialTotalEnergy*EnergyRatio))*
    ((1-EnergyRatio)/(2.0*initialTotalEnergy*EnergyRatio));

  // Calculate rejection function at 0, 1 and Emax
  G4double gfunction0 = RejectionFunction(0);
  G4double gfunction1 = RejectionFunction(1);
  G4double gfunctionEmax = RejectionFunction(gMaxEnergy);


  // Calculate Maximum value 
  G4double gMaximum = std::max(gfunction0,gfunction1);
  gMaximum = std::max(gMaximum,gfunctionEmax);

  G4double rand, gfunctionTest, randTest;

  do{
    rand = G4UniformRand();
    rand = rand/(1-rand+1.0/gMaxEnergy);
    gfunctionTest = RejectionFunction(rand);
    randTest = G4UniformRand();

  }while(randTest > (gfunctionTest/gMaximum));

  theta = sqrt(rand)/initialTotalEnergy;


  return theta;
}
//

G4double G4Generator2BS::RejectionFunction(G4double value) const
{

  G4double argument = (1+value)*(1+value);

  G4double gfunction = (4+log(rejection_argument3+(z/argument)))*
    ((4*EnergyRatio*value/argument)-rejection_argument1)+rejection_argument2;

  return gfunction;

}

void G4Generator2BS::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Bremsstrahlung Angular Generator is 2BS Generator from 2BS Koch & Motz distribution (Rev Mod Phys 31(4), 920 (1959))" << G4endl;
  G4cout << "Sampling algorithm adapted from PIRS-0203" << G4endl;
  G4cout << "\n" << G4endl;
} 

