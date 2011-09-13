 
#ifndef G4LENDFission_h
#define G4LENDFission_h 1

// Class Description
// Manager of LEND (Low Energy Nuclear Data) target (nucleus) 
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

#include "G4LENDModel.hh"

class G4LENDFission : public G4LENDModel
{

   public: 
  
     G4LENDFission( G4ParticleDefinition* pd )
     :G4LENDModel( "LENDFission" ) 
     { 
        proj = pd; 
       
//        theModelName = "LENDFission for "; 
//        theModelName += proj->GetParticleName(); 
        create_used_target_map();
     };
  
     ~G4LENDFission(){;};
  
     G4HadFinalState* ApplyYourself( const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus );
};

#endif
