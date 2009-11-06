// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: 
//	 	1 July 2009 creation by L. Desorgher based on a modification of G4He3 		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint particles are used in Reverse/Adjoint Monte Carlo simulations. New adjoint 
//		processes act on adjoint particles when they are  tracked backward in the geometry. 
//		The use of adjoint particles instead of "normal" particles during a reverse simulation 
//		is based on an idea of M. Asai.   
//
#ifndef G4AdjointHe3_h
#define G4AdjointHe3_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4AdjointIons.hh"

// ######################################################################
// ###                            ADJOINT He3                         ###
// ######################################################################

class G4AdjointHe3 : public G4AdjointIons
{
 private:
   static G4AdjointHe3* theInstance;
   G4AdjointHe3(){}
   ~G4AdjointHe3(){}

 public:
   static G4AdjointHe3* Definition();
   static G4AdjointHe3* He3Definition();
   static G4AdjointHe3* He3();
};

#endif

