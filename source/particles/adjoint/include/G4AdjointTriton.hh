// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: 
//	 	1 July 2009 creation by L. Desorgher based on a modification of G4Triton		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint particles are used in Reverse/Adjoint Monte Carlo simulations. New adjoint 
//		processes act on adjoint particles when they are  tracked backward in the geometry. 
//		The use of adjoint particles instead of "normal" particles during a reverse simulation 
//		is based on an idea of M. Asai.   
//
#ifndef G4AdjointTriton_h
#define G4AdjointTriton_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4AdjointIons.hh"

// ######################################################################
// ###                          ADJOINT TRITON                        ###
// ######################################################################

class G4AdjointTriton : public G4AdjointIons
{
 private:
   static G4AdjointTriton* theInstance;
   G4AdjointTriton(){}
   ~G4AdjointTriton(){}

 public:
   static G4AdjointTriton* Definition();
   static G4AdjointTriton* TritonDefinition();
   static G4AdjointTriton* Triton();
};

#endif
