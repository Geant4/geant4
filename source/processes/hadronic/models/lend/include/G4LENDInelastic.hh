#ifndef G4LENDInelastic_h
#define G4LENDInelastic_h 1

// Class Description
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

class G4LENDInelastic : public G4LENDModel
{

   public: 
  
     G4LENDInelastic( G4ParticleDefinition* pd )
     { 
        proj = pd; 
       
//        theModelName = "LENDInelastic for "; 
//        theModelName += proj->GetParticleName(); 
        create_used_target_map();
     };
  
     ~G4LENDInelastic(){;};
  
     G4HadFinalState* ApplyYourself( const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus );
};

#endif
