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
/// \file eventgenerator/pythia/pythia8decayer/src/SingleParticleGun.cc
/// \brief Implementation of the SingleParticleGun class
///
/// \author J. Yarba; FNAL

#include "SingleParticleGun.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

// SPECIAL CASE - needed for treatment of polarization for tau's
#include "G4TauMinus.hh"
#include "G4TauPlus.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SingleParticleGun::SingleParticleGun( const G4String& pname, const double pmom )
   : G4VUserPrimaryGeneratorAction(),
     fGun(nullptr),
     fMomentum(pmom)
{

   int nparts = 1;
   fGun = new G4ParticleGun( nparts );

   G4ParticleTable* pdt = G4ParticleTable::GetParticleTable();
   G4ParticleDefinition* pd = pdt->FindParticle( pname );
   fGun->SetParticleDefinition( pd );
   fGun->SetParticlePosition( G4ThreeVector(0.,0.,0.) );
   double mass = pd->GetPDGMass();
   double energy = std::sqrt( fMomentum*fMomentum + mass*mass );
   fGun->SetParticleEnergy( energy );

   // NOTE: do not set momentum directions 
   //       as it'll be generated event-by-event

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SingleParticleGun::~SingleParticleGun()
{

   delete fGun;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SingleParticleGun::GeneratePrimaries( G4Event* evt )
{

   if ( !fGun )
   {
      G4cout << " ParticleGun is NOT set up; BAIL OUT from this event " << G4endl;
      G4cout << " Check if you are using ctor " 
             << "SinglePartGun(const G4string& partname, const double pmom ) " 
             << G4endl;
      return;
   }

   double phi   = -1. * CLHEP::pi + CLHEP::twopi * G4UniformRand();
   double theta = CLHEP::pi/4. + CLHEP::halfpi * G4UniformRand();

   double x = std::sin(theta) * std::cos(phi);
   double y = std::sin(theta) * std::sin(phi);
   double z = std::cos(theta);

   fGun->SetParticleMomentumDirection( G4ThreeVector(x,y,z) );

   // SPECIAL CASE: for particle such as tau polarization should be -1 * p3
   //               while for tau_bar it should be +1
   //
   // NOTE: one could probably compare by name but that'd mean comparing strings...
   //
   if ( fGun->GetParticleDefinition() == G4TauMinus::TauMinus() )
   {
      fGun->SetParticlePolarization( G4ThreeVector(-1.0*x,-1.0*y,-1.0*z) );
   }
   else if ( fGun->GetParticleDefinition() == G4TauPlus::TauPlus() )
   {
      fGun->SetParticlePolarization( G4ThreeVector(x,y,z) );
   }

   fGun->GeneratePrimaryVertex(evt);

   return;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
