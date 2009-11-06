// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: 
//	 	1 July 2009 creation by L. Desorgher based on a modification of G4Ions		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint particles are used in Reverse/Adjoint Monte Carlo simulations. New adjoint 
//		processes act on adjoint particles when they are  tracked backward in the geometry. 
//		The use of adjoint particles instead of "normal" particles during a reverse simulation 
//		is based on an idea of M. Asai.   
//

#ifndef G4AdjointIons_h
#define G4AdjointIons_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

// ######################################################################
// ###                     ADJOINT Ions                               ###
// ######################################################################

class G4AdjointIons : public G4ParticleDefinition
{
 // Class Description
 //  This is the base class for all nuclei including pre-defined 
 //  light nuclei such as deuteron, alpha, and proton (Hydrogen) 
 //  All nuclei/ions created on the fly are objects of this class
 //  Atomic number and atomic mass are vaild only for particles derived
 //  from this class.  This class has Excitation Energy in addition to
 //  the normal particle properties.

 protected:
   G4AdjointIons(){};


 public: //With Description
   G4AdjointIons(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable,  G4bool              shortlived,
       const G4String&     subType ="",
       G4int               anti_encoding =0,
       G4double            excitation = 0.0
   );

 public:
   virtual    			~G4AdjointIons();
   G4AdjointIons*    			IonsDefinition();
   G4AdjointIons*    			Ions();

 public:  //With Description
   // Get excitation energy of nucleus
   G4double GetExcitationEnergy() const ; 
  
  private:
   G4double theExcitationEnergy; 

};

inline
 G4AdjointIons* G4AdjointIons::Ions() 
{
  return this;
}

inline
 G4double G4AdjointIons::GetExcitationEnergy() const 
{
  return theExcitationEnergy;
}

#endif








