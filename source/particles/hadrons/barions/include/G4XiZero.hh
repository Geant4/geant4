// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4XiZero.hh,v 1.3 1999-10-03 09:11:52 kurasige Exp $
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

// Each class inheriting from G4VBaryon
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4XiZero_h
#define G4XiZero_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VBaryon.hh"

// ######################################################################
// ###                          XiZero                                ###
// ######################################################################

class G4XiZero : public G4VBaryon
{
 private:
   static G4XiZero theXiZero;
   static G4double  theXiZeroLengthCut;
   static G4double* theXiZeroKineticEnergyCuts;

 private:
   G4XiZero(
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
   virtual ~G4XiZero(){}

   static G4XiZero* XiZeroDefinition();
   static G4XiZero* XiZero() {return &theXiZero;}
   static G4double GetCuts() {return theXiZeroLengthCut;}   
   static G4double* GetCutsInEnergy() {return theXiZeroKineticEnergyCuts;};

   virtual void SetCuts(G4double aCut); 
};

#endif
