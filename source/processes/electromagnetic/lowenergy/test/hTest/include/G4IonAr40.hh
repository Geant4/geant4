#ifndef G4IonAr40_h
#define G4IonAr40_h 1

//---------------------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//---------------------------------------------------------------------------
//
// ClassName:   G4IonAr40
//
// Class Description:
// The new static ion Ar40+ is defined as G4VIon.
// Each class inheriting from G4VIon
// corresponds to a particle type; one and only one
// instance for each class is guaranteed.
// Class Description - end
//  
// Authors:   08.04.01  V.Ivanchenko 
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4ios.hh"
#include "G4VIon.hh"

// ######################################################################
// ###                          IonAr40                                 ###
// ######################################################################

class G4IonAr40 : public G4VIon
{
 private:
   static G4IonAr40 theIonAr40;
   static G4double  theIonAr40LengthCut;
   static G4double* theIonAr40KineticEnergyCuts;

public: // Without description
   G4IonAr40(
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

   virtual ~G4IonAr40();
  
   static G4IonAr40*    IonAr40Definition();
   static G4IonAr40*    IonAr40() {return &theIonAr40;}
   static G4double GetCuts() {return theIonAr40LengthCut;}   
   static G4double* GetCutsInEnergy() {return theIonAr40KineticEnergyCuts;};

   void SetCuts(G4double aCut); 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4IonAr40::SetCuts(G4double aCut)
{
  CalcEnergyCuts(aCut);
  theIonAr40LengthCut = theCutInMaxInteractionLength;  
  theIonAr40KineticEnergyCuts = theKineticEnergyCuts;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
