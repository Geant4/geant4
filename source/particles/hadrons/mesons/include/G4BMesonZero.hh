// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BMesonZero.hh,v 1.2 1999-04-13 08:37:11 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//
//      Created,             Hisaya Kurashige, 15 June 1997
// **********************************************************************
//  Change both methods to get the pointer into non-inlined H.Kurashige 4 Aug. 1998
// ----------------------------------------------------------------

// Each class inheriting from G4VMeson
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4BMesonZero_h
#define G4BMesonZero_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VMeson.hh"

// ######################################################################
// ###                         BMesonZero                             ###
// ######################################################################

class G4BMesonZero : public G4VMeson
{
 private:
   static G4BMesonZero theBMesonZero;
   static G4double  theBMesonZeroLengthCut;
   static G4double* theBMesonZeroKineticEnergyCuts;

 private: // constructors are hide as private  
   G4BMesonZero(
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
   virtual ~G4BMesonZero(){}

   static G4BMesonZero* BMesonZeroDefinition();
   static G4BMesonZero* BMesonZero();
   static G4double GetCuts() {return theBMesonZeroLengthCut;}   
   static G4double* GetCutsInEnergy() {return theBMesonZeroKineticEnergyCuts;};

   virtual void SetCuts(G4double aCut); 
};

#endif
