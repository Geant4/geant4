#ifndef G4LENDCrossSection_h
#define G4LENDCrossSection_h 1

// Class Description
// Cross Sections for a LEND (Low Energy Nuclear Data)
// LEND is Geant4 interface for GIDI (General Interaction Data Interface) 
// which gives a discription of nuclear and atomic reactions, such as
//    Binary collision cross sections
//    Particle number multiplicity distributions of reaction products
//    Energy and angular distributions of reaction products
//    Derived calculational constants
// GIDI is developped at Lawrence Livermore National Laboratory
// Class Description - End

// 071025 First implementation done by T. Koi (SLAC/SCCS)
// 101118 Name modifications for release T. Koi (SLAC/PPA)

#include "G4LENDHeader.hh"

#include "G4LENDManager.hh"
#include "G4LENDUsedTarget.hh"

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsTable.hh"

#include <map>

class G4LENDCrossSection : public G4VCrossSectionDataSet
{
   public:
   
      G4LENDCrossSection(const G4String name = "" );
   
      ~G4LENDCrossSection();
   
      G4bool IsApplicable( const G4DynamicParticle* , const G4Element* );

      void BuildPhysicsTable( const G4ParticleDefinition& );

      void DumpPhysicsTable( const G4ParticleDefinition& );

      G4double GetCrossSection( const G4DynamicParticle* , const G4Element* , G4double aT );

      G4double GetIsoCrossSection( const G4DynamicParticle* , const G4Isotope* , G4double );

      G4double GetZandACrossSection(const G4DynamicParticle* , G4int /*Z*/, G4int /*A*/, G4double aTemperature);

      void ChangeDefaultEvaluation( G4String name ){ default_evaluation = name; recreate_used_target_map(); };
      //void AllowNaturalAbundanceTarget(){ allow_nat = true; recreate_used_target_map(); };
      //void AllowAnyCandidateTarget(){ allow_any = true; recreate_used_target_map(); };
      void AllowNaturalAbundanceTarget(){ allow_nat = true; };
      void AllowAnyCandidateTarget(){ allow_any = true; };

      //Hadronic Framework still does not handle isotope in GPIL
      //G4VDiscreteProcess::PostStepGetPhysicalInteractionLenght()
      // G4HadronicProcess::GetMeanFreePath()
      //  G4double GetCrossSection(const G4DynamicParticle*, const G4Material*)
      //   G4double GetCrossSection(const G4DynamicParticle*, const G4Element*, G4double aTemperature); 
      G4bool IsIsoApplicable( const G4DynamicParticle* , G4int /*ZZ*/, G4int /*AA*/) { return true; }

   private:
   
      virtual G4double getLENDCrossSection( G4GIDI_target* , G4double , G4double ) { return 0.0; };

      std::map< G4int , G4LENDUsedTarget* > usedTarget_map;

      G4String default_evaluation;
      G4bool allow_nat;
      G4bool allow_any;

      G4LENDManager* lend_manager;
      void recreate_used_target_map();

   protected :
      //G4String name;
      G4ParticleDefinition* proj;
      void create_used_target_map();

};
#endif
