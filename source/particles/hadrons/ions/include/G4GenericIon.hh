// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GenericIon.hh,v 1.4 2001-03-12 05:45:43 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      4-th Dec 1998, H.Kurashige
// ****************************************************************
// This class is used only by G4IonTable and not for tracking
// G4IonTable creates various ions other than alpha,deuteron,triton,
// and He3. Processes for these ions will be same as ones for 
// this "GenericIon". So, user should register processes for ions
// to this class in his/her UserPhysicsList

#ifndef G4GenericIon_h
#define G4GenericIon_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VIon.hh"

// ######################################################################
// ###                          GenericIon                            ###
// ######################################################################

class G4GenericIon : public G4VIon
{
 private:
   static G4GenericIon theGenericIon;
   static G4double  theGenericIonLengthCut;
   static G4double* theGenericIonKineticEnergyCuts;

 public:
   G4GenericIon(
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
   virtual ~G4GenericIon();

   static G4GenericIon*    GenericIonDefinition();
   static G4GenericIon*    GenericIon(){return &theGenericIon;}
   static G4double GetCuts() {return theGenericIonLengthCut;}   
   static G4double* GetCutsInEnergy() {return theGenericIonKineticEnergyCuts;};

   void SetCuts(G4double aCut); 
   virtual void RestoreCuts(G4double cutInLength,
			    const G4double* cutInEnergy );
};

inline void G4GenericIon::SetCuts(G4double aCut)
{
  CalcEnergyCuts(aCut);
  theGenericIonLengthCut = theCutInMaxInteractionLength;  
  theGenericIonKineticEnergyCuts = theKineticEnergyCuts;
  
}

inline void G4GenericIon::RestoreCuts(G4double cutInLength,
			    const G4double* cutInEnergy )
{
  G4ParticleWithCuts::RestoreCuts(cutInLength, cutInEnergy);
  theGenericIonLengthCut = theCutInMaxInteractionLength;  
  theGenericIonKineticEnergyCuts = theKineticEnergyCuts;
  
}
#endif
