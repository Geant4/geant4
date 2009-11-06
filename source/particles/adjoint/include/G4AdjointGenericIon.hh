// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: 
//	 	10 July 2009 creation by L. Desorgher based on a modification of G4GenericIon		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint particles are used in Reverse/Adjoint Monte Carlo simulations. New adjoint 
//		processes act on adjoint particles when they are  tracked backward in the geometry. 
//		The use of adjoint particles instead of "normal" particles during a reverse simulation 
//		is based on an idea of M. Asai.   
//

#ifndef G4AdjointGenericIon_h
#define G4AdjointGenericIon_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4AdjointIons.hh"

// ######################################################################
// ###                          GenericIon                            ###
// ######################################################################

class G4AdjointGenericIon : public G4AdjointIons
{
 private:
   static G4AdjointGenericIon* theInstance;
   G4AdjointGenericIon(){}
   ~G4AdjointGenericIon(){}

 public:
   static G4AdjointGenericIon* Definition();
   static G4AdjointGenericIon* GenericIonDefinition();
   static G4AdjointGenericIon* GenericIon();
};

#endif
