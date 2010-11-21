#ifndef G4LENDFissionCrossSection_h
#define G4LENDFissionCrossSection_h 1

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

#include "G4LENDCrossSection.hh"

class G4LENDFissionCrossSection : public G4LENDCrossSection
{

   public:
   
      G4LENDFissionCrossSection()
      {;};
      G4LENDFissionCrossSection( G4ParticleDefinition* pd )
      {
         proj = pd; 
         name = "LENDFission for ";
         name += proj->GetParticleName();
         create_used_target_map();
      };
   
      ~G4LENDFissionCrossSection(){;};
   
   private:
//                                                       ke         temperature
      G4double getLENDCrossSection( GIDI4GEANT_target* , G4double , G4double );

};
#endif
