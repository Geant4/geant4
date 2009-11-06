// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: 
//	 	1 July 2009 creation by L. Desorgher based on a modification of G4Deuteron 		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint particles are used in Reverse/Adjoint Monte Carlo simulations. New adjoint 
//		processes act on adjoint particles when they are  tracked backward in the geometry. 
//		The use of adjoint particles instead of "normal" particles during a reverse simulation 
//		is based on an idea of M. Asai.   
//
#ifndef G4AdjointDeuteron_h
#define G4AdjointDeuteron_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4AdjointIons.hh"

// ######################################################################
// ###                          ADJOINT DEUTERON                      ###
// ######################################################################

class G4AdjointDeuteron : public G4AdjointIons
{
 private:
   static G4AdjointDeuteron* theInstance;
   G4AdjointDeuteron(){}
   ~G4AdjointDeuteron(){}

 public:
   static G4AdjointDeuteron* Definition();
   static G4AdjointDeuteron* DeuteronDefinition();
   static G4AdjointDeuteron* Deuteron();
};

#endif






