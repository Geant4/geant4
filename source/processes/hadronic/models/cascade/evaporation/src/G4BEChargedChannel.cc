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
// Implementation of the HETC88 code into Geant4.
// Evaporation and De-excitation parts
// T. Lampen, Helsinki Institute of Physics, May-2000
//
// 20120608  M. Kelsey -- Change vars 's','m','m2' to avoid name collisions

#include "G4BEChargedChannel.hh"
#include "G4SystemOfUnits.hh"

G4BEChargedChannel::G4BEChargedChannel()
{
    verboseLevel = 0;
}


G4BEChargedChannel::~G4BEChargedChannel()
{
}


void G4BEChargedChannel::calculateProbability()
{
  G4int residualZ = nucleusZ - particleZ;
  G4int residualA = nucleusA - particleA;

//    Check if nucleus is too small, if this evaporation channel
//    leads to an impossible residual nucleus or if there is no enough
//    energy.
  if ( nucleusA < 2.0 * particleA || 
       nucleusZ < 2.0 * particleZ ||
       residualA <= residualZ || 
       excitationEnergy - getThresh() - correction < 0 )
    {
      if ( verboseLevel >= 6 )
	G4cout << "G4BEChargedChannel : calculateProbability for " << getName() << " = 0 " << G4endl;
      emissionProbability = 0;
      return;
    }

  // In HETC88 s-s0 was used in std::exp( s ), in which s0 was either 50 or
  // max(s_i), where i goes over all channels.

  G4double levelParam = getLevelDensityParameter();
  G4double slevel = 2 * std::sqrt( levelParam  * ( excitationEnergy - getThresh() - correction ) );
  G4double constant = A / 2 * ( 2 * spin + 1 ) * ( 1 + coulombFactor() );
  G4double eye1     = ( slevel*slevel - 3 * slevel + 3 ) / ( 4 * levelParam*levelParam ) * std::exp( slevel );

  emissionProbability = constant * std::pow( G4double(residualA), 0.6666666 ) * eye1;

  if ( verboseLevel >= 6 )
    G4cout << "G4BEChargedChannel : calculateProbability for " << getName() << G4endl
	   << "                    res A = " << residualA << G4endl 
	   << "                    res Z = " << residualZ << G4endl 
	   << "                 c*(c_i+1) = "<< constant << G4endl
	   << "                  qmfactor = "<< qmFactor() << G4endl
	   << "             coulombfactor = "<< coulombFactor() << G4endl
	   << "                        E = " << excitationEnergy << G4endl
	   << "               correction = " << correction << G4endl
	   << "                     eye1 = " << eye1 << G4endl
	   << "               levelParam = " << levelParam << G4endl
	   << "                   thresh = " << getThresh() << G4endl
	   << "                        s = " << s << G4endl
	   << "              probability = " << emissionProbability << G4endl;

  return;
}


G4double G4BEChargedChannel::sampleKineticEnergy()
{
  G4double levelParam;
  levelParam = getLevelDensityParameter();
  
  const G4double xMax  = excitationEnergy - getThresh() - correction; // maximum number
  const G4double xProb = ( - 1 + std::sqrt ( 1 + 4 * levelParam * xMax ) ) / ( 2 * levelParam ); // most probable value
  const G4double maxProb = xProb * std::exp ( 2 * std::sqrt ( levelParam * ( xMax - xProb ) ) ); // maximum value of P(x)

  // Sample x according to density function P(x) with rejection method
  G4double r1;
  G4double r2;
  G4int koe=0;
  do
    {
      r1 = G4UniformRand() * xMax;
      r2 = G4UniformRand() * maxProb;
      koe++;
    }
  while (  r1 * std::exp ( 2 * std::sqrt ( levelParam * ( xMax - r1 ) ) )  < r2 );

//  G4cout << "Q ch " << koe << G4endl;
  G4double kineticEnergy = r1 + getCoulomb(); // add coulomb potential;

  if ( verboseLevel >= 10 )
    G4cout << " G4BENeutronChannel : sampleKineticEnergy() " << G4endl
	   << "       kinetic n e = " << kineticEnergy << G4endl
	   << "        levelParam = " << levelParam << G4endl
	   << "             thresh= " << getThresh() << G4endl;

  return kineticEnergy;
}


G4double G4BEChargedChannel::coulombFactorForProton()
{
  //  Coefficient c_p:s for empirical cross section formula are
  //  defined with the proton constant.  See Dostrovsky, Phys. Rev.,
  //  vol. 116, 1959.
  G4double t[7] = { 0.08 , 0 , -0.06 , -0.1 , -0.1 , -0.1 , -0.1 };
  G4int Z = nucleusZ - particleZ;

  if ( Z >= 70.0 ) return t[6];
  if ( Z <= 10.0 ) return t[0];
  
  // Linear interpolation
  G4int   n = G4int( 0.1 * Z + 1.0 ); 
  G4float x = ( 10 * n - Z ) * 0.1; 
  G4double ret_val =  x * t[n - 2] + ( 1.0 - x ) * t[n-1];

  return ret_val;
}


G4double G4BEChargedChannel::qmFactorForProton()
{
  //  Coefficient k_p for empirical cross section formula are defined
  //  with the proton constant.  See Dostrovsky, Phys. Rev., vol. 116,
  //  1959
  G4double t[7] = {  0.36, 0.51, 0.60, 0.66, 0.68, 0.69, 0.69 };
  G4int Z = nucleusZ - particleZ;

  if ( Z >= 70.0 ) return t[6];
  if ( Z <= 10.0 ) return t[0];

  // Linear interpolation
  G4int   n = G4int( 0.1 * Z + 1.0 ); 
  G4float x = ( 10 * n - Z ) * 0.1; 
  return x * t[n - 2] + ( 1.0 - x ) * t[n-1];
}


G4double G4BEChargedChannel::qmFactorForAlpha()
{
//  Coefficient k_alpha for empirical cross section formula presented
//  in Dostrovsky, Phys. Rev., vol. 116, 1959

  G4double t[7] = {  0.77, 0.81, 0.85, 0.89, 0.93, 0.97,  1.00 };
  G4int Z = nucleusZ - particleZ;

  if ( Z >= 70.0 ) return t[6];
  if ( Z <= 10.0 ) return t[0];

  // Linear interpolation
  G4int   n = G4int( 0.1 * Z + 1.0 ); 
  G4float x = ( 10 * n - Z ) * 0.1; 
  return x * t[n - 2] + ( 1.0 - x ) * t[n-1];
}
