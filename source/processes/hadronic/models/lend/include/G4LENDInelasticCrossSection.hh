#ifndef G4LENDInelasticCrossSection_h
#define G4LENDInelasticCrossSection_h 1

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

#include "G4LENDCrossSection.hh"

class G4LENDInelasticCrossSection : public G4LENDCrossSection
{

   public:
   
      G4LENDInelasticCrossSection()
      {;};
      G4LENDInelasticCrossSection( G4ParticleDefinition* pd )
      {
         proj = pd; 
         name = "LENDInelastic for ";
         name += proj->GetParticleName();
         create_used_target_map();
      };
   
      ~G4LENDInelasticCrossSection(){;};
   
   private:
//                                                        ke         temperature
      G4double getLENDCrossSection( GIDI4GEANT_target* , G4double , G4double );

};
#endif
