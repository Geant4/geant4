// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EtaPrime.hh,v 1.3 1999-12-15 14:51:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, 8 June 1998 Hisaya Kurashige
// **********************************************************************
// ------------------------------------------------------------

// Each class inheriting from G4VMeson
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4EtaPrime_h
#define G4EtaPrime_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VMeson.hh"

// ######################################################################
// ###                         EtaPrime                               ###
// ######################################################################

class G4EtaPrime : public G4VMeson
{
 private:
   static G4EtaPrime  theEtaPrime;
   static G4double    theEtaPrimeLengthCut;
   static G4double*   theEtaPrimeKineticEnergyCuts;

 private: // constructors are hide as private  
   G4EtaPrime(
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
   virtual  ~G4EtaPrime(){}

   static G4EtaPrime*      EtaPrimeDefinition();
   static G4EtaPrime*      EtaPrime(){return &theEtaPrime;}
   static G4double GetCuts() {return theEtaPrimeLengthCut;}   
   static G4double* GetCutsInEnergy() {return theEtaPrimeKineticEnergyCuts;};

   virtual void        SetCuts(G4double aCut);
};

#endif





