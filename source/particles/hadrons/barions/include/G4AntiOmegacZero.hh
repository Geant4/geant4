// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AntiOmegacZero.hh,v 1.2 1999-04-13 08:31:49 kurasige Exp $
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

#ifndef G4AntiOmegacZero_h
#define G4AntiOmegacZero_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VBarion.hh"

// ######################################################################
// ###                          AntiOmegacZero                            ###
// ######################################################################

class G4AntiOmegacZero : public G4VBarion
{
 private:
   static G4AntiOmegacZero theAntiOmegacZero;
   static G4double  theAntiOmegacZeroLengthCut;
   static G4double* theAntiOmegacZeroKineticEnergyCuts;

 private:
   G4AntiOmegacZero(
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
   virtual ~G4AntiOmegacZero(){}

   static G4AntiOmegacZero* AntiOmegacZeroDefinition();
   static G4AntiOmegacZero* AntiOmegacZero();
   static G4double GetCuts() {return theAntiOmegacZeroLengthCut;}   
   static G4double* GetCutsInEnergy() {return theAntiOmegacZeroKineticEnergyCuts;};

   virtual void SetCuts(G4double aCut); 
};


#endif
