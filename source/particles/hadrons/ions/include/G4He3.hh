// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4He3.hh,v 1.2 1999-04-13 08:24:06 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      24-th April 1998,H.Kurashige
// ****************************************************************
// ----------------------------------------------------------------------
// Each class inheriting from G4VIon
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4He3_h
#define G4He3_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VIon.hh"

// ######################################################################
// ###                          He3                                 ###
// ######################################################################

class G4He3 : public G4VIon
{
 private:
   static G4He3 theHe3;
   static G4double  theHe3LengthCut;
   static G4double* theHe3KineticEnergyCuts;

 public:
   G4He3(
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
   virtual ~G4He3();

   static G4He3*    He3Definition();
   static G4He3*    He3() {return &theHe3;}
   static G4double GetCuts() {return theHe3LengthCut;}   
   static G4double* GetCutsInEnergy() {return theHe3KineticEnergyCuts;};

   void SetCuts(G4double aCut); 
};

inline void G4He3::SetCuts(G4double aCut)
{
  CalcEnergyCuts(aCut);
  theHe3LengthCut = theCutInMaxInteractionLength;  
  theHe3KineticEnergyCuts = theKineticEnergyCuts;
  
}

#endif

