// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SigmaPlus.hh,v 1.4 1999-12-15 14:50:55 gunter Exp $
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

#ifndef G4SigmaPlus_h
#define G4SigmaPlus_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VBaryon.hh"

// ######################################################################
// ###                          SigmaPlus                             ###
// ######################################################################

class G4SigmaPlus : public G4VBaryon
{
 private:
   static G4SigmaPlus theSigmaPlus;
   static G4double  theSigmaPlusLengthCut;
   static G4double* theSigmaPlusKineticEnergyCuts;

 private:
   G4SigmaPlus(
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
   virtual ~G4SigmaPlus(){}

   static G4SigmaPlus* SigmaPlusDefinition();
   static G4SigmaPlus* SigmaPlus() {return &theSigmaPlus;}
   static G4double GetCuts() {return theSigmaPlusLengthCut;}   
   static G4double* GetCutsInEnergy() {return theSigmaPlusKineticEnergyCuts;};

   virtual void SetCuts(G4double aCut); 
};

inline void G4SigmaPlus::SetCuts(G4double aCut)
{
  CalcEnergyCuts(aCut);
  theSigmaPlusLengthCut = theCutInMaxInteractionLength;  
  theSigmaPlusKineticEnergyCuts = theKineticEnergyCuts;
  
}


#endif
