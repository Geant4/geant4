// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Positron.hh,v 1.4 2001-03-12 05:47:38 kurasige Exp $
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
//  Added SetCuts, L.Urban, 12 June 1996
//  Revised, Hisaya Kurashige, 4 July 1996
//  Revised, Hisaya Kurashige, 7 July 1996
//  Added not static GetEnergyCuts() and GetLengthCuts(), G.Cosmo, 11 July 1996
//  Revised, Hisaya Kurashige, 15 Dec 1996
// ----------------------------------------------------------------
// Each class inheriting from G4VLepton
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.

#ifndef G4Positron_h
#define G4Positron_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VLepton.hh"
#include "G4Electron.hh"

// ######################################################################
// ###                         POSITRON                               ###
// ######################################################################

class G4Positron : public G4VLepton
{
 private:
   static G4Positron thePositron;
   static G4double  thePositronLengthCut;
   static G4double* thePositronKineticEnergyCuts;
        
 protected:  
   G4double ComputeLoss(G4double AtomicNumber, G4double KineticEnergy) const;

   void BuildRangeVector( const G4Material* aMaterial,
			  const G4LossTable* aLossTable,
			  G4double       maxEnergy,     
			  G4double       aMass,
                          G4RangeVector* rangeVector
                         );

   friend void G4Electron::BuildRangeVector(
				  const G4Material* aMaterial,
				  const G4LossTable* aLossTable,
				  G4double       maxEnergy,     
				  G4double       aMass,
                                  G4RangeVector* rangeVector
                         );


 private:   //hide constructor as private 
   G4Positron(	
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
   virtual ~G4Positron(){};	
 
   static G4Positron* PositronDefinition();
   static G4Positron* Positron();
   static G4double GetCuts() {return thePositronLengthCut;}
   static G4double* GetCutsInEnergy() {return thePositronKineticEnergyCuts;};

   virtual void SetCuts(G4double aCut); 
   virtual void RestoreCuts(G4double cutInLength,
			    const G4double* cutInEnergy );
};

inline void G4Positron::SetCuts(G4double aCut)
{
  CalcEnergyCuts(aCut);
  thePositronLengthCut = theCutInMaxInteractionLength;  
  thePositronKineticEnergyCuts = theKineticEnergyCuts;
}

inline void G4Positron::RestoreCuts(G4double cutInLength,
			    const G4double* cutInEnergy )
{
  G4ParticleWithCuts::RestoreCuts(cutInLength, cutInEnergy);
  thePositronLengthCut = theCutInMaxInteractionLength;  
  thePositronKineticEnergyCuts = theKineticEnergyCuts;
}

inline G4Positron*  G4Positron::Positron()
{  return &thePositron; }

#endif



