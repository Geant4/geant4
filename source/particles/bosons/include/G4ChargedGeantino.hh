// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ChargedGeantino.hh,v 1.1 1999-01-07 16:09:47 gunter Exp $
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
//  Added particle definitions, H.Kurashige, 19 April 1996
//  Revised, G.Cosmo, 6 June 1996
//  Added not static GetEnergyCuts() and GetLengthCuts(), G.Cosmo, 11 July 1996
// ----------------------------------------------------------------

// Each class inheriting from G4VBoson
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4ChargedGeantino_h
#define G4ChargedGeantino_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VBoson.hh"

// ######################################################################
// ###                         ChargedGeantino                        ###
// ######################################################################

class G4ChargedGeantino : public G4VBoson
{
 private:
   static G4ChargedGeantino theChargedGeantino;
   static G4double  theChargedGeantinoLengthCut;
   static G4double* theChargedGeantinoKineticEnergyCuts;

 private:
   G4ChargedGeantino(
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
   static G4ChargedGeantino* ChargedGeantinoDefinition();
   static G4ChargedGeantino* ChargedGeantino();
   static G4double  GetCuts() {return theChargedGeantinoLengthCut;}   
   static G4double* GetCutsInEnergy() {return theChargedGeantinoKineticEnergyCuts;};

   void   SetCuts(G4double aCut);
};

inline G4ChargedGeantino* G4ChargedGeantino::ChargedGeantino()
{  return &theChargedGeantino; }

#endif
