// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4KaonMinus.hh,v 1.4 2001-03-12 05:45:45 kurasige Exp $
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

// Each class inheriting from G4VMeson
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4KaonMinus_h
#define G4KaonMinus_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VMeson.hh"

// ######################################################################
// ###                        KAONMINUS                               ###
// ######################################################################

class G4KaonMinus : public G4VMeson
{
 private:
   static G4KaonMinus theKaonMinus;
   static G4double  theKaonMinusLengthCut;
   static G4double* theKaonMinusKineticEnergyCuts;

 private: // constructors are hide as private  
   G4KaonMinus(
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
   virtual ~G4KaonMinus(){}

   static G4KaonMinus* KaonMinusDefinition();
   static G4KaonMinus* KaonMinus() {return &theKaonMinus;}
   static G4double GetCuts() {return theKaonMinusLengthCut;}   
   static G4double* GetCutsInEnergy() {return theKaonMinusKineticEnergyCuts;};

   virtual void SetCuts(G4double aCut); 
   virtual void RestoreCuts(G4double cutInLength,
			    const G4double* cutInEnergy );
};

inline void G4KaonMinus::SetCuts(G4double aCut)
{
  G4ParticleWithCuts::SetCuts(aCut);
  theKaonMinusLengthCut = theCutInMaxInteractionLength;  
  theKaonMinusKineticEnergyCuts = theKineticEnergyCuts;
}

inline void G4KaonMinus::RestoreCuts(G4double cutInLength,
			    const G4double* cutInEnergy )
{
  G4ParticleWithCuts::RestoreCuts(cutInLength, cutInEnergy);
  theKaonMinusLengthCut = theCutInMaxInteractionLength;  
  theKaonMinusKineticEnergyCuts = theKineticEnergyCuts;
}

#endif
