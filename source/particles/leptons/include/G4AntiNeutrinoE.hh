// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AntiNeutrinoE.hh,v 1.2 1999-04-13 08:20:16 kurasige Exp $
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

// Each class inheriting from G4VLepton
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4AntiNeutrinoE_h
#define G4AntiNeutrinoE_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VLepton.hh"

// ######################################################################
// ###                       ANTI NEUTRINO                            ###
// ######################################################################

class G4AntiNeutrinoE : public G4VLepton
{
 private:
   static G4AntiNeutrinoE theAntiNeutrinoE;
   static G4double  theAntiNeutrinoELengthCut;
   static G4double* theAntiNeutrinoEKineticEnergyCuts;

 private: // constructors are hide as private  
   G4AntiNeutrinoE(
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
   virtual ~G4AntiNeutrinoE(){}
 
   static G4AntiNeutrinoE* AntiNeutrinoEDefinition();
   static G4AntiNeutrinoE* AntiNeutrinoE(){return &theAntiNeutrinoE;}
   static G4double  GetCuts() {return theAntiNeutrinoELengthCut;}   
   static G4double* GetCutsInEnergy() {return theAntiNeutrinoEKineticEnergyCuts;};

   virtual void SetCuts(G4double aCut); 
};

#endif

