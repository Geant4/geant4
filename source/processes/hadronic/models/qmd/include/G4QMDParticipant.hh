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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:    G4QMDParticipant.hh 
//
//      Author: Koi, Tatsumi (tkoi@slac.stanford.edu)       
// 
//      Creation date: 29 March 2007
// -----------------------------------------------------------------------------
//
// 081120 Add hit flag and related methods

#ifndef G4QMDParticipant_hh
#define G4QMDParticipant_hh

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"

class G4QMDParticipant 
{
   public:
                                                // momentum      position
      G4QMDParticipant( const G4ParticleDefinition* , G4ThreeVector , G4ThreeVector );
      ~G4QMDParticipant();

      void SetDefinition( const G4ParticleDefinition* pd ) { definition = pd; };
      const G4ParticleDefinition* GetDefinition() { return definition; };

      void SetPosition( G4ThreeVector r ) { position = r; };
      G4ThreeVector GetPosition() { return position; };

      void SetMomentum( G4ThreeVector p ) { momentum = p; };
      G4ThreeVector GetMomentum() { return momentum; };

      G4double GetMass() { return definition->GetPDGMass()/CLHEP::GeV; };

      G4LorentzVector Get4Momentum();

      G4double GetKineticEnergy() { return Get4Momentum().e() - GetMass(); };

      G4int GetBaryonNumber() { return definition->GetBaryonNumber(); };
      G4int GetNuc() { return definition->GetBaryonNumber(); };

      G4int GetChargeInUnitOfEplus() { return int ( definition->GetPDGCharge()/CLHEP::eplus ); };

      void UnsetInitialMark() { projectile = false; target = false; }
      void UnsetHitMark() { hit = false; }
      G4bool IsThisHit() { return hit; }
      void SetHitMark() { hit = true; }

      void SetProjectile() { projectile = true; }
      void SetTarget() { target = true; }
      G4bool IsThisProjectile() { return projectile; }
      G4bool IsThisTarget() { return target; }

   private:
      const G4ParticleDefinition* definition;
      G4ThreeVector momentum;
      G4ThreeVector position;

      G4bool projectile; 
      G4bool target; 
      G4bool hit; 
};

#endif
