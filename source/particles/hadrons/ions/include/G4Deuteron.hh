// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Deuteron.hh,v 1.1 1999-01-07 16:10:05 gunter Exp $
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

// Each class inheriting from G4VIon
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4Deuteron_h
#define G4Deuteron_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VIon.hh"

// ######################################################################
// ###                          DEUTERON                              ###
// ######################################################################

class G4Deuteron : public G4VIon
{
 private:
   static G4Deuteron theDeuteron;
   static G4double  theDeuteronLengthCut;
   static G4double* theDeuteronKineticEnergyCuts;

 public:
   G4Deuteron(
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
   static G4Deuteron* DeuteronDefinition();
   static G4Deuteron* Deuteron(){return &theDeuteron;}
   static G4double GetCuts() {return theDeuteronLengthCut;}   
   static G4double* GetCutsInEnergy() {return theDeuteronKineticEnergyCuts;};

   void SetCuts(G4double aCut); 
};

inline void G4Deuteron::SetCuts(G4double aCut)
{
  CalcEnergyCuts(aCut);
  theDeuteronLengthCut = theCutInMaxInteractionLength;  
  theDeuteronKineticEnergyCuts = theKineticEnergyCuts;
  
}

#endif
