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

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4Nucleus.hh"  
#include "G4BEGammaDeexcitation.hh"


G4BEGammaDeexcitation::G4BEGammaDeexcitation()
{
}


G4BEGammaDeexcitation::~G4BEGammaDeexcitation()
{
}


void G4BEGammaDeexcitation::setVerboseLevel( G4int level )
{
  verboseLevel = level ;
}


void G4BEGammaDeexcitation::setNucleusA( G4int a )
{
  nucleusA = a;
}


void G4BEGammaDeexcitation::setNucleusZ( G4int z )
{
  nucleusZ = z;
}


void G4BEGammaDeexcitation::setExcitationEnergy( G4double energy )
{
  excitationEnergy = energy;
}


G4double  G4BEGammaDeexcitation::sampleKineticEnergy()
{
  return excitationEnergy  * G4UniformRand(); 
}


G4DynamicParticle *  G4BEGammaDeexcitation::emit()
{
  // Isotropic distribution assumed to gammas
  G4double u, v, w;
  G4DynamicParticle * pParticle = new G4DynamicParticle;
  pParticle -> SetDefinition( G4Gamma::Gamma() );
  pParticle -> SetKineticEnergy( sampleKineticEnergy() ); 
  isotropicCosines( u, v, w );
  pParticle -> SetMomentumDirection( u, v, w );  
  return pParticle;
}


void G4BEGammaDeexcitation::isotropicCosines( G4double & u, 
					      G4double & v, 
					      G4double & w )
{
  // Samples isotropic random direction cosines.
  G4double CosTheta = 1.0 - 2.0 * G4UniformRand();
  G4double SinTheta = std::sqrt( 1.0 - CosTheta * CosTheta );
  G4double Phi = twopi * G4UniformRand();

  u = std::cos( Phi ) * SinTheta;
  v = std::cos( Phi ) * CosTheta,
  w = std::sin( Phi );

  return;
}
