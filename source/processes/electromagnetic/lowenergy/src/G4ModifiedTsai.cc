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
// File name:     G4ModifiedTsai
//
// Author:        Andreia Trindade (andreia@lip.pt)
//                Pedro Rodrigues  (psilva@lip.pt)
//                Luis Peralta     (luis@lip.pt)
// 
// Creation date: 21 March 2003
//
// Modifications: 
// 21 Mar 2003       A. Trindade    First implementation acording with new design
// 24 Mar 2003                      & Fix in Tsai generator in order to prevent theta generation above pi
//
// Class Description: 
//
// Concrete base class for Bremsstrahlung Angular Distribution Generation - Tsai Model
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//    

#include "G4ModifiedTsai.hh"
#include "Randomize.hh"
//    

G4ModifiedTsai::G4ModifiedTsai(const G4String& name):G4VBremAngularDistribution(name)
{;}

//    

G4ModifiedTsai::~G4ModifiedTsai() 
{;}

//

G4double G4ModifiedTsai::PolarAngle(const G4double initial_energy,
				    const G4double, // final_energy
				    const G4int ) // Z
{

  // Sample gamma angle (Z - axis along the parent particle).
  // Universal distribution suggested by L. Urban (Geant3 manual (1993) 
  // Phys211) derived from Tsai distribution (Rev Mod Phys 49,421(1977))

  G4double totalEnergy = initial_energy + electron_mass_c2;   

  const G4double a1 = 0.625, a2 = 3.*a1, d = 27.;
  G4double u, theta = 0;

  do{
  u = - std::log(G4UniformRand()*G4UniformRand());

  if (9./(9.+d) > G4UniformRand()) u /= a1;
  else                             u /= a2;
    
  theta = u*electron_mass_c2/totalEnergy;
  }while(u > (totalEnergy*pi/electron_mass_c2));

  return theta;
}
//

void G4ModifiedTsai::PrintGeneratorInformation() const
{

  G4cout << "\n" << G4endl;
  G4cout << "Bremsstrahlung Angular Generator is Modified Tsai" << G4endl;
  G4cout << "Universal distribution suggested by L. Urban (Geant3 manual (1993) Phys211)" << G4endl;
  G4cout << "Derived from Tsai distribution (Rev Mod Phys 49,421(1977)) \n" << G4endl;
} 
