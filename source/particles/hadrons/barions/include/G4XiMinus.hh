// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4XiMinus.hh,v 1.1 1999-01-07 16:09:56 gunter Exp $
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
// ----------------------------------------------------------------

// Each class inheriting from G4VBarion
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4XiMinus_h
#define G4XiMinus_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VBarion.hh"

// ######################################################################
// ###                          XiMinus                               ###
// ######################################################################

class G4XiMinus : public G4VBarion
{
 private:
   static G4XiMinus theXiMinus;
   static G4double  theXiMinusLengthCut;
   static G4double* theXiMinusKineticEnergyCuts;

 private:
   G4XiMinus(
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
   static G4XiMinus* XiMinusDefinition();
   static G4XiMinus* XiMinus() {return &theXiMinus;}
   static G4double GetCuts() {return theXiMinusLengthCut;}   
   static G4double* GetCutsInEnergy() {return theXiMinusKineticEnergyCuts;};

   void SetCuts(G4double aCut); 
};

inline void G4XiMinus::SetCuts(G4double aCut)
{
  CalcEnergyCuts(aCut);
  theXiMinusLengthCut = theCutInMaxInteractionLength;  
  theXiMinusKineticEnergyCuts = theKineticEnergyCuts;
  
}


#endif
