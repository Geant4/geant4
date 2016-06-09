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
// $Id: G4DipBustGenerator.cc,v 1.1 2010-10-14 15:17:48 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4DipBustGenerator
//
// Author:  Vladimir Grichine     
// 
// Creation date: 17 May 2011
//
// Modifications: 
// 
//
// Class Description: 
//
// Bremsstrahlung Angular Distribution Generation 
// suggested  the dipole approximation in the rest frame of electron 
// busted in the laboratory frame.
//
// Class Description: End 
//
// -------------------------------------------------------------------
//

#include "G4DipBustGenerator.hh"
#include "Randomize.hh"

G4DipBustGenerator::G4DipBustGenerator(const G4String&)
  : G4VBremAngularDistribution("AngularGenUrban")
{}    

G4DipBustGenerator::~G4DipBustGenerator() 
{}

G4double G4DipBustGenerator::PolarAngle(const G4double eTkin,
				    const G4double, // final_energy
				    const G4int ) // Z
{
  G4double c, cosTheta, delta, cofA, signc = 1., a, power = 1./3.;
  G4double gamma, beta, theta;

  c = 4. - 8.*G4UniformRand();
  a = c;
 
  if( c < 0. )
  {
    signc = -1.;
    a     = -c;
  }
  delta  = std::sqrt(a*a+4.);
  delta += a;
  delta *= 0.5; 

  cofA = -signc*std::pow(delta, power);

  cosTheta = cofA - 1./cofA;

  gamma = 1. + eTkin/electron_mass_c2;
  beta = std::sqrt(1. - 1./gamma/gamma);

  cosTheta = (cosTheta + beta)/(1 + cosTheta*beta);

  theta = std::acos(cosTheta);

  if( theta < 0. )  theta = 0.;
  if( theta > pi )  theta = pi;
  // G4cout <<"theta = "<<theta<<"; ";

  return theta;
}

void G4DipBustGenerator::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Bremsstrahlung Angular Generator is Modified Tsai" << G4endl;
  G4cout << "Distribution suggested by L.Urban (Geant3 manual (1993) Phys211)" 
	 << G4endl;
  G4cout << "Derived from Tsai distribution (Rev Mod Phys 49,421(1977)) \n" 
	 << G4endl;
} 
