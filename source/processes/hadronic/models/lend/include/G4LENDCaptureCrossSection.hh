#ifndef G4LENDCaptureCrossSection_h
#define G4LENDCaptureCrossSection_h 1

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

class G4LENDCaptureCrossSection : public G4LENDCrossSection
{

   public:
   
      //G4LENDCaptureCrossSection()
      //{;};
      G4LENDCaptureCrossSection( G4ParticleDefinition* pd )
      :G4LENDCrossSection("LENDCaptureCrossSection")
      {
         proj = pd; 
         //name = "LEND Capture Cross Section for ";
         //name += proj->GetParticleName();
         //create_used_target_map();
      };
   
      ~G4LENDCaptureCrossSection(){;};
   
   private:
//                                                        ke         temperature
      G4double getLENDCrossSection( G4GIDI_target* , G4double , G4double );

};
#endif
