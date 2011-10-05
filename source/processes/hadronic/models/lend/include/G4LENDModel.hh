#ifndef G4LENDModel_h
#define G4LENDModel_h 1

// Class Description
// Final state production model for a LEND (Low Energy Nuclear Data) 
// LEND is Geant4 interface for GIDI (General Interaction Data Interface) 
// which gives a discription of nuclear and atomic reactions, such as
//    Binary collision cross sections
//    Particle number multiplicity distributions of reaction products
//    Energy and angular distributions of reaction products
//    Derived calculational constants
// LEND is developped at Lawrence Livermore National Laboratory
// Class Description - End

// 071025 First implementation done by T. Koi (SLAC/SCCS)
// 101118 Name modifications for release T. Koi (SLAC/PPA)

#include "G4LENDHeader.hh"
#include "G4LENDManager.hh"
#include "G4LENDUsedTarget.hh"
#include "G4HadronicInteraction.hh"
#include "globals.hh"

class G4LENDModel : public G4HadronicInteraction
{

   public: 
  
      //G4LENDModel();
      G4LENDModel( G4String name="LENDModel" );
      ~G4LENDModel();
  
      virtual G4HadFinalState * ApplyYourself( const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus );

      void ChangeDefaultEvaluation( G4String name ){ default_evaluation = name; recreate_used_target_map(); };
      void AllowNaturalAbundanceTarget(){ allow_nat = true; recreate_used_target_map(); };
      void AllowAnyCandidateTarget(){ allow_any = true; recreate_used_target_map(); };
      //Same argument to the CrossSectionDataSet 
      void BuildPhysicsTable( const G4ParticleDefinition& ){ recreate_used_target_map(); };

   private:

      G4String default_evaluation;
      G4bool allow_nat;
      G4bool allow_any;

   protected:

      void create_used_target_map();
      void recreate_used_target_map();

      G4ParticleDefinition* proj;
      G4LENDManager* lend_manager;
      std::map< G4int , G4LENDUsedTarget* > usedTarget_map;
};

#endif
