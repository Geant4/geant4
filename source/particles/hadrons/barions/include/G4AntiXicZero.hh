// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AntiXicZero.hh,v 1.4 1999-12-15 14:50:55 gunter Exp $
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

// Each class inheriting from G4VBaryon
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4AntiXicZero_h
#define G4AntiXicZero_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VBaryon.hh"

// ######################################################################
// ###                          AntiXicZero                                ###
// ######################################################################

class G4AntiXicZero : public G4VBaryon
{
 private:
   static G4AntiXicZero theAntiXicZero;
   static G4double  theAntiXicZeroLengthCut;
   static G4double* theAntiXicZeroKineticEnergyCuts;

 private:
   G4AntiXicZero(
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
   virtual ~G4AntiXicZero(){}
 
   static G4AntiXicZero* AntiXicZeroDefinition();
   static G4AntiXicZero* AntiXicZero();
   static G4double GetCuts() {return theAntiXicZeroLengthCut;}   
   static G4double* GetCutsInEnergy() {return theAntiXicZeroKineticEnergyCuts;};

   virtual void SetCuts(G4double aCut); 
};

#endif
