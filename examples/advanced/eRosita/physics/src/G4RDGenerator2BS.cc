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
// GEANT4 Class file
//
//
// File name:     G4RDGenerator2BS
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

#include "G4RDGenerator2BS.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
//    

G4RDGenerator2BS::G4RDGenerator2BS(const G4String& name):G4RDVBremAngularDistribution(name)
{;}

//    

G4RDGenerator2BS::~G4RDGenerator2BS() 
{;}

//

G4double G4RDGenerator2BS::PolarAngle(const G4double initial_energy,
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
  G4double gMaxEnergy = (pi*initialTotalEnergy)*(pi*initialTotalEnergy);

  G4double Zeff = std::sqrt(static_cast<G4double>(Z) * (static_cast<G4double>(Z) + 1.0));
  z = (0.00008116224*(std::pow(Zeff,0.3333333)));

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

  theta = std::sqrt(rand)/initialTotalEnergy;


  return theta;
}
//

G4double G4RDGenerator2BS::RejectionFunction(G4double value) const
{

  G4double argument = (1+value)*(1+value);

  G4double gfunction = (4+std::log(rejection_argument3+(z/argument)))*
    ((4*EnergyRatio*value/argument)-rejection_argument1)+rejection_argument2;

  return gfunction;

}

void G4RDGenerator2BS::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Bremsstrahlung Angular Generator is 2BS Generator from 2BS Koch & Motz distribution (Rev Mod Phys 31(4), 920 (1959))" << G4endl;
  G4cout << "Sampling algorithm adapted from PIRS-0203" << G4endl;
  G4cout << "\n" << G4endl;
} 

