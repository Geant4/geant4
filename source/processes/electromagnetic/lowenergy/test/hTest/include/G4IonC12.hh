#ifndef G4IonC12_h
#define G4IonC12_h 1

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
// ClassName:   G4IonC12
//
// Class Description:
// The new static ion C12+ is defined as G4VIon.
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4IonC12::SetCuts(G4double aCut)
{
  CalcEnergyCuts(aCut);
  theIonC12LengthCut = theCutInMaxInteractionLength;  
  theIonC12KineticEnergyCuts = theKineticEnergyCuts;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
