// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PionMinus.hh,v 1.4 2001-03-12 05:45:45 kurasige Exp $
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

#ifndef G4PionMinus_h
#define G4PionMinus_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VMeson.hh"

// ######################################################################
// ###                        PIONMINUS                               ###
// ######################################################################

class G4PionMinus : public G4VMeson
{
 private:
   static G4PionMinus thePionMinus;
   static G4double  thePionMinusLengthCut;
   static G4double* thePionMinusKineticEnergyCuts;

 private: // constructors are hide as private  
   G4PionMinus(
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
  virtual      ~G4PionMinus(){}
   static      G4PionMinus* PionMinusDefinition();
   static      G4PionMinus* PionMinus(){return &thePionMinus;}
   static G4double GetCuts() {return thePionMinusLengthCut;}   
   static G4double* GetCutsInEnergy() {return thePionMinusKineticEnergyCuts;};

   virtual void SetCuts(G4double aCut); 
   virtual void RestoreCuts(G4double cutInLength,
			    const G4double* cutInEnergy );
};

inline void G4PionMinus::SetCuts(G4double aCut)
{
  G4ParticleWithCuts::SetCuts(aCut);
  thePionMinusLengthCut = theCutInMaxInteractionLength;  
  thePionMinusKineticEnergyCuts = theKineticEnergyCuts;
}

inline void G4PionMinus::RestoreCuts(G4double cutInLength,
			    const G4double* cutInEnergy )
{
  G4ParticleWithCuts::RestoreCuts(cutInLength, cutInEnergy);
  thePionMinusLengthCut = theCutInMaxInteractionLength;  
  thePionMinusKineticEnergyCuts = theKineticEnergyCuts;
}

#endif
