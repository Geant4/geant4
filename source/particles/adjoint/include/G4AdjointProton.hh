
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: 
//	 	13 February 2009 creation by L. Desorgher based on a modification of G4Proton 		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint particles are used in Reverse/Adjoint Monte Carlo simulations. New adjoint 
//		processes act on adjoint particles when they are  tracked backward in the geometry. 
//		The use of adjoint particles instead of "normal" particles during a reverse simulation 
//		is based on an idea of M. Asai.   
//
#ifndef G4AdjointProton_h
#define G4AdjointProton_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4AdjointIons.hh"

// ######################################################################
// ###                          ADJOINT PROTON                        ###
// ######################################################################

class G4AdjointProton : public G4ParticleDefinition
{
 private:
   static G4AdjointProton* theInstance;
   G4AdjointProton(){}
   ~G4AdjointProton(){}

 public:
   static G4AdjointProton* Definition();
   static G4AdjointProton* AdjointProtonDefinition();
   static G4AdjointProton* AdjointProton();
};


#endif
