// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OmegacZero.hh,v 1.2 1999-04-13 08:32:57 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      4-th April 1996, G.Cosmo
// ****************************************************************
//  Added particle definitions, H.Kurashige, 14 Feb 19
//  Change both methods to get the pointer into non-inlined H.Kurashige 4 Aug. 1998
// ----------------------------------------------------------------

// Each class inheriting from G4VBarion
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4OmegacZero_h
#define G4OmegacZero_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VBarion.hh"

// ######################################################################
// ###                          OmegacZero                            ###
// ######################################################################

class G4OmegacZero : public G4VBarion
{
 private:
   static G4OmegacZero theOmegacZero;
   static G4double  theOmegacZeroLengthCut;
   static G4double* theOmegacZeroKineticEnergyCuts;

 private:
   G4OmegacZero(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable
   );

 public:
   virtual  ~G4OmegacZero(){}

   static G4OmegacZero* OmegacZeroDefinition();
   static G4OmegacZero* OmegacZero();
   static G4double GetCuts() {return theOmegacZeroLengthCut;}   
   static G4double* GetCutsInEnergy() {return theOmegacZeroKineticEnergyCuts;};

   virtual void SetCuts(G4double aCut); 
};


#endif
