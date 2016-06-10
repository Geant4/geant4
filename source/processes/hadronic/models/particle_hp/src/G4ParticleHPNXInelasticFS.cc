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
// particle_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPNXInelasticFS.hh"
#include "G4Nucleus.hh"
#include "G4Alpha.hh"

G4HadFinalState * G4ParticleHPNXInelasticFS::ApplyYourself(const G4HadProjectile & theTrack)
{
// these are the particle types in the final state

  G4ParticleDefinition * theDefs[2];
  theDefs[0] = G4Neutron::Neutron();
  
// fill the final state  
  G4ParticleHPInelasticBaseFS::BaseApply(theTrack, theDefs, 1);
  
// return the result
   return theResult.Get();
}

void G4ParticleHPNXInelasticFS::
Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String & aFSType, G4ParticleDefinition* projectile)
{
   G4ParticleHPInelasticBaseFS::Init(A, Z, M, dirName, aFSType, projectile);
   G4double ResidualA = 0;
   G4double ResidualZ = 0;
   if( projectile == G4Neutron::Neutron() ) {
     ResidualA = A-5;
     ResidualZ = Z-2;
   } else if( projectile == G4Proton::Proton() ) {
     ResidualA = A-5;
     ResidualZ = Z-1;
   } else if( projectile == G4Deuteron::Deuteron() ) {
     ResidualA = A-4;
     ResidualZ = Z-1;
   } else if( projectile == G4Triton::Triton() ) {
     ResidualA = A-3;
     ResidualZ = Z-1;
   } else if( projectile == G4He3::He3() ) {
     ResidualA = A-3;
     ResidualZ = Z;
   } else if( projectile == G4Alpha::Alpha() ) {
     ResidualA = A-2;
     ResidualZ = Z;
   }

   //   G4double ResidualA = A-5; // wrong, to be improves @@@
   //   G4double ResidualZ = Z-2; // dito @@@
   G4ParticleHPInelasticBaseFS::InitGammas(ResidualA, ResidualZ);
}
