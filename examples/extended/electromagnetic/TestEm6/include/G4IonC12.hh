// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IonC12.hh,v 1.2 1999-10-28 02:15:14 vnivanch Exp $
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

// Class Description:
// The new ion C12+ is defined as G4VIon.
// Each class inheriting from G4VIon
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.
// Class Description - end


#ifndef G4IonC12_h
#define G4IonC12_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VIon.hh"

// ######################################################################
// ###                          IonC12                                 ###
// ######################################################################

class G4IonC12 : public G4VIon
{
 private:
   static G4IonC12  theIonC12;
   static G4double  theIonC12LengthCut;
   static G4double* theIonC12KineticEnergyCuts;

public: // Without description

   G4IonC12(
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
   virtual ~G4IonC12();
  
   static G4IonC12*    IonC12Definition();
   static G4IonC12*    IonC12() {return &theIonC12;}
   static G4double GetCuts() {return theIonC12LengthCut;}   
   static G4double* GetCutsInEnergy() {return theIonC12KineticEnergyCuts;};

   void SetCuts(G4double aCut); 
};

inline void G4IonC12::SetCuts(G4double aCut)
{
  CalcEnergyCuts(aCut);
  theIonC12LengthCut = theCutInMaxInteractionLength;  
  theIonC12KineticEnergyCuts = theKineticEnergyCuts;
  
}

#endif
