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
// Implementation of the HETC88 code into Geant4.
// Evaporation and De-excitation parts
// T. Lampen, Helsinki Institute of Physics, May-2000
//
// 20120608  M. Kelsey -- Change vars 's','m','m2' to avoid name collisions

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4ParticleTable.hh"    
#include "G4Nucleus.hh"  
#include "G4BENeutronChannel.hh"


G4BENeutronChannel::G4BENeutronChannel()
{
  name         = "neutron";
  particleA    = 1;
  particleZ    = 0;
  verboseLevel = 0;
  rho          = 0;
}


G4BENeutronChannel::~G4BENeutronChannel()
{
}


void G4BENeutronChannel::calculateProbability()
{
  const G4int residualZ = nucleusZ - particleZ;
  const G4int residualA = nucleusA - particleA;

  if ( nucleusA  <  2.0 * particleA ||  
       nucleusZ  <  2.0 * particleZ ||  
       residualA <= residualZ       ||  
       excitationEnergy - getThresh() - correction < 0 )
    {
      if ( verboseLevel >= 6 )
	G4cout << "G4BENeutronChannel : calculateProbability = 0 " << G4endl;
      emissionProbability = 0;
      return;
    }

  // In HETC88 s-s0 was used in std::exp( s ), in which s0 was either 50 or
  // max(s_i), where i goes over all channels.

  const G4double levelParam = getLevelDensityParameter();
  
  const G4double slevel = 2 * std::sqrt( levelParam  * ( excitationEnergy - getThresh() - correction ) );
  const G4double eye0 = std::exp( slevel ) * ( slevel - 1 ) /  ( 2 * levelParam );
  const G4double eye1 = ( slevel*slevel - 3*slevel +3 ) * std::exp( slevel ) / ( 4 * levelParam*levelParam ) ;
  
  emissionProbability = std::pow( G4double(residualA), 0.666666 ) * alpha() * ( eye1 + beta() * eye0 );
  
  if ( verboseLevel >= 6 )
    G4cout << "G4BENeutronChannel : calculateProbability " << G4endl
	   << "                    res A = " << residualA << G4endl 
	   << "                    res Z = " << residualZ << G4endl 
	   << "                    alpha = " << alpha() << G4endl 
	   << "                     beta = " << beta() << G4endl
	   << "                        E = " << excitationEnergy << G4endl
	   << "               correction = " << correction << G4endl
	   << "                     eye1 = " << eye1 << G4endl
	   << "                     eye0 = " <<  eye0 << G4endl
	   << "               levelParam = " << levelParam << G4endl
	   << "                   thresh = " << getThresh() << G4endl
	   << "                        s = " << s << G4endl
	   << "              probability = " << emissionProbability << G4endl;

  return;
}


G4double  G4BENeutronChannel::sampleKineticEnergy()
{
  ////////////////
  // A random number is sampled from the density function
  // P(x) = x * std::exp ( 2 std::sqrt ( a ( xMax - x ) ) )  [not normalized],
  // x belongs to [ 0, xMax ]
  // with the 'Hit or Miss' -method
  // Kinetic energy is this energy scaled properly
  
  G4double levelParam;
  levelParam = getLevelDensityParameter();
  
  const G4double xMax  = excitationEnergy - getThresh() - correction + beta(); // maximum number
  const G4double xProb = ( - 1 + std::sqrt ( 1 + 4 * levelParam * xMax ) ) / ( 2 * levelParam ); // most probable value
  const G4double maxProb = xProb * std::exp ( 2 * std::sqrt ( levelParam * ( xMax - xProb ) ) ); // maximum value of P(x)

  // Sample x according to density function P(x) with rejection method
  G4double r1;
  G4double r2;
  G4int koe=0;
  do
    {
      r1 = beta() + G4UniformRand() * ( xMax - beta() );
      r2 = G4UniformRand() * maxProb;
      koe++;
    }
  while (  r1 * std::exp ( 2 * std::sqrt ( levelParam * ( xMax - r1 ) ) )  < r2 );

//  G4cout <<  koe << G4endl;
  G4double kineticEnergy = r1 - beta();

  if ( verboseLevel >= 10 )
    G4cout << " G4BENeutronChannel : sampleKineticEnergy() " << G4endl
	   << "       kinetic n e = " << kineticEnergy << G4endl
	   << "        levelParam = " << levelParam << G4endl
	   << "             thresh= " << getThresh() << G4endl
	   << "               beta= " << beta() << G4endl;

  return kineticEnergy;
}


G4DynamicParticle *  G4BENeutronChannel::emit()
{
  G4double u;
  G4double v;
  G4double w;
  G4DynamicParticle * pParticle = new G4DynamicParticle;

  pParticle -> SetDefinition( G4Neutron::Neutron() );
  pParticle -> SetKineticEnergy( sampleKineticEnergy() ); 
  isotropicCosines( u, v, w );
  pParticle -> SetMomentumDirection( u , v , w );  

  return pParticle;
}


G4double G4BENeutronChannel::alpha()
{
  const G4double residualA = nucleusA - particleA;
  return 0.76 + 1.93 * std::pow( residualA, -0.33333 );
}


G4double G4BENeutronChannel::beta()
{
  G4double residualA = nucleusA - particleA;
  return ( 1.66 * std::pow ( residualA, -0.66666 ) - 0.05 )/alpha()*MeV;
}

