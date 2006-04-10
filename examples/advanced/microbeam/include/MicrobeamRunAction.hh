// -------------------------------------------------------------------
// $Id: MicrobeamRunAction.hh,v 1.2 2006-04-10 14:47:31 sincerti Exp $
// -------------------------------------------------------------------

#ifndef MicrobeamRunAction_h
#define MicrobeamRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4ThreeVector.hh"

#include "MicrobeamDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

class MicrobeamRunAction : public G4UserRunAction
{
public:
  
  MicrobeamRunAction(MicrobeamDetectorConstruction*);
  ~MicrobeamRunAction();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
    
  void  SetRndmFreq(G4int   val)  {saveRndm = val;}
  G4int GetRndmFreq()             {return saveRndm;}

  void AddDoseN(G4float dose){ DoseN += dose;}
  void SetDoseN(G4float dose){ DoseN = dose;}
  G4float GetDoseN(){return DoseN;}

  void AddDoseC(G4float dose){ DoseC += dose;}
  void SetDoseC(G4float dose){ DoseC = dose;}
  G4float GetDoseC(){return DoseC;}

  G4int GetNumEvent(){return numEvent;}
  void SetNumEvent(G4int i){numEvent = i;}

  G4int GetNbOfHitsGas(){return nbOfHitsGas;}
  void AddNbOfHitsGas(){nbOfHitsGas = nbOfHitsGas+1;}

  void SetMassPhantom(G4float mP){ massPhantom = mP;}
  G4float GetMassPhantom(){return massPhantom;}

  void SetMassNucleus(G4float mN){ massNucleus = mN;}
  G4float GetMassNucleus(){return massNucleus;}

  void SetMassCytoplasm(G4float mC){ massCytoplasm = mC;}
  G4float GetMassCytoplasm(){return massCytoplasm;}

  void AddDoseBox(G4int i, G4float x){ dose3DDose[i] +=x;}
  G4float GetDoseBox(G4int i){ return dose3DDose[i];}
  
  G4ThreeVector GetVectCell(G4int i) {return mapVoxels[i];}

private:

  MicrobeamDetectorConstruction* Detector;    
  MicrobeamPhantomConfiguration myMicrobeamPhantomConfiguration;  

  G4int saveRndm;
  G4int numEvent;
  G4int nbOfPixels;
  G4int nbOfHitsGas;
  G4float SP;
  G4float R;
  G4float RnElec;
  G4float RcElec;
  G4float DoseN;
  G4float DoseC;
  G4float DoseNElec;
  G4float DoseCElec;
  G4float massPhantom;
  G4float massCytoplasm;
  G4float massNucleus;
  G4bool boolSP;
  
  G4int * x3DDose;
  G4int * y3DDose;
  G4int * z3DDose;
  G4float * dose3DDose;
  G4ThreeVector * mapVoxels;

};

#endif
