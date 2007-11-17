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

#ifndef G4QMDParticipant_hh
#define G4QMDParticipant_hh

#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"

class G4QMDParticipant 
{
   public:
                                                // momentum      position
      G4QMDParticipant( G4ParticleDefinition* , G4ThreeVector , G4ThreeVector );
      ~G4QMDParticipant();

      void SetDefinition( G4ParticleDefinition* pd ) { definition = pd; };
      G4ParticleDefinition* GetDefinition() { return definition; };

      void SetPosition( G4ThreeVector r ) { position = r; };
      G4ThreeVector GetPosition() { return position; };

      void SetMomentum( G4ThreeVector p ) { momentum = p; };
      G4ThreeVector GetMomentum() { return momentum; };

      G4double GetMass() { return definition->GetPDGMass()/GeV; };

      G4LorentzVector Get4Momentum();

      G4double GetKineticEnergy() { return Get4Momentum().e() - GetMass(); };

      G4int GetBaryonNumber() { return definition->GetBaryonNumber(); };
      G4int GetNuc() { return definition->GetBaryonNumber(); };

      G4int GetChargeInUnitOfEplus() { return int ( definition->GetPDGCharge()/eplus ); };

      void UnsetInitialMark() { projectile = false; target = false; }
      void SetProjectile() { projectile = true; }
      void SetTarget() { target = true; }
      G4bool IsThisProjectile() { return projectile; }
      G4bool IsThisTarget() { return target; }

   private:
      G4ParticleDefinition* definition;
      G4ThreeVector momentum;
      G4ThreeVector position;

      G4bool projectile; 
      G4bool target; 
};

#endif
