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

#include "G4BEDeuteronChannel.hh"
#include "G4Deuteron.hh"


G4BEDeuteronChannel::G4BEDeuteronChannel()
{
  name = "deuteron";
  verboseLevel = 0;
  particleA = 2;
  particleZ = 1;
  rho = 0.70588;
  exmass = G4NucleiProperties::GetMassExcess( particleA, particleZ ); 
  A=2;
  spin = 1;
}


G4BEDeuteronChannel::~G4BEDeuteronChannel()
{
}


G4DynamicParticle*  G4BEDeuteronChannel::emit()
{
  G4double u, v, w;
  G4DynamicParticle * pParticle = new G4DynamicParticle;
  pParticle -> SetDefinition( G4Deuteron::Deuteron() );
  pParticle -> SetKineticEnergy( sampleKineticEnergy() ); 
  isotropicCosines( u, v, w );
  pParticle -> SetMomentumDirection( u, v, w );  
  return pParticle;
}


G4double G4BEDeuteronChannel::coulombFactor()
{
  // Coefficient c_i representing the variation of charged-particle
  // capture cross sections with Z for each particle.  See Dostrovsky,
  // Phys. Rev. 116, 1959.
  return 0.5 * coulombFactorForProton();
}


G4double G4BEDeuteronChannel::qmFactor()
{
  // Coefficient k_i representing the quantum mechanical barrier
  // penetration, see Dostrovsky, Phys. Rev. 116, 1959.
  return qmFactorForProton() + 0.06;
}
