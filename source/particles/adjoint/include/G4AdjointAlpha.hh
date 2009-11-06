// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: 
//	 	1 July 2009 creation by L. Desorgher based on a modification of G4Alpha 		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint particles are used in Reverse/Adjoint Monte Carlo simulations. New adjoint 
//		processes act on adjoint particles when they are  tracked backward in the geometry. 
//		The use of adjoint particles instead of "normal" particles during a reverse simulation 
//		is based on an idea of M. Asai.   
//
#ifndef G4AdjointAlpha_h
#define G4AdjointAlpha_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4AdjointIons.hh"

// ######################################################################
// ###                          ADJOINT ALPHA                         ###
// ######################################################################

class G4AdjointAlpha : public G4AdjointIons
{
 private:
   static G4AdjointAlpha* theInstance;
   G4AdjointAlpha(){}
   ~G4AdjointAlpha(){}

 public:
   static G4AdjointAlpha* Definition();
   static G4AdjointAlpha* AlphaDefinition();
   static G4AdjointAlpha* Alpha();
};

#endif
