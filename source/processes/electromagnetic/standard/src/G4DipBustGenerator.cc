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
// $Id: G4DipBustGenerator.cc 74581 2013-10-15 12:03:25Z gcosmo $
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
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

G4DipBustGenerator::G4DipBustGenerator(const G4String&)
  : G4VEmAngularDistribution("DipBustGen")
{}    

G4DipBustGenerator::~G4DipBustGenerator() 
{}

G4ThreeVector& 
G4DipBustGenerator::SampleDirection(const G4DynamicParticle* dp,
				    G4double, G4int, const G4Material*)
{
  G4double a, c, cosTheta, delta, cofA, signc = 1.;

  G4double eTkin = dp->GetKineticEnergy();

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

  cofA = -signc*G4Exp(G4Log(delta)/3.0);

  cosTheta = cofA - 1./cofA;

  G4double tau = eTkin/electron_mass_c2;
  G4double beta = std::sqrt(tau*(tau + 2.))/(tau + 1.);

  cosTheta = (cosTheta + beta)/(1 + cosTheta*beta);

  G4double sinTheta = std::sqrt((1 - cosTheta)*(1 + cosTheta));
  G4double phi  = twopi*G4UniformRand(); 

  fLocalDirection.set(sinTheta*std::cos(phi), sinTheta*std::sin(phi),cosTheta);
  fLocalDirection.rotateUz(dp->GetMomentumDirection());

  return fLocalDirection;

}

G4double G4DipBustGenerator::PolarAngle(const G4double eTkin,
				    const G4double, // final_energy
				    const G4int ) // Z
{
  G4double c, cosTheta, delta, cofA, signc = 1., a;
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

  cofA = -signc*G4Exp(G4Log(delta)/3.0);

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
  G4cout << "Angular Generator based on classical formula from" << G4endl;
  G4cout << "J.D. Jackson, Classical Electrodynamics, Wiley, New York 1975" 
	 << G4endl;
} 
