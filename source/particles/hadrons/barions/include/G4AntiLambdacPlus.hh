// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AntiLambdacPlus.hh,v 1.3 1999-10-03 09:11:21 kurasige Exp $
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

#ifndef G4AntiLambdacPlus_h
#define G4AntiLambdacPlus_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VBaryon.hh"

// ######################################################################
// ###                          AntiLambdacPlus                       ###
// ######################################################################

class G4AntiLambdacPlus : public G4VBaryon
{
 private:
   static G4AntiLambdacPlus theAntiLambdacPlus;
   static G4double  theAntiLambdacPlusLengthCut;
   static G4double* theAntiLambdacPlusKineticEnergyCuts;

 private:
   G4AntiLambdacPlus(
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
   virtual ~G4AntiLambdacPlus(){}

   static G4AntiLambdacPlus* AntiLambdacPlusDefinition();
   static G4AntiLambdacPlus* AntiLambdacPlus();
   static G4double GetCuts() {return theAntiLambdacPlusLengthCut;}   
   static G4double* GetCutsInEnergy() {return theAntiLambdacPlusKineticEnergyCuts;};

   virtual void SetCuts(G4double aCut); 
};

inline void G4AntiLambdacPlus::SetCuts(G4double aCut)
{
  CalcEnergyCuts(aCut);
  theAntiLambdacPlusLengthCut = theCutInMaxInteractionLength;  
  theAntiLambdacPlusKineticEnergyCuts = theKineticEnergyCuts;
  
}
#endif
