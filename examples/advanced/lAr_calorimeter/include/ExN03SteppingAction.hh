// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN03SteppingAction.hh,v 1.1 2002-10-01 13:43:56 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef ExN03SteppingAction_h
#define ExN03SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ExN03SteppingAction : public G4UserSteppingAction
{
  public:
    ExN03SteppingAction();
    virtual ~ExN03SteppingAction();

    virtual void UserSteppingAction(const G4Step*);

  private:
  G4int EventNo;
  G4int IDnow, IDout, IDold;
  G4int NTracks, NSecondaries;
  G4double OutOfWorldTracksData[6000][11], Secondaries[6000][11];
  G4double EdepFCALEm, EdepFCALHad;

  G4ThreeVector PrimaryVertex, PrimaryDirection;
  G4ThreeVector SecondaryVertex, SecondaryDirection;
  G4ThreeVector Distance, VectorProduct, VectorProductNorm;
  G4double VectorProductMagnitude, DistOfClosestApproach;

  public:
  void initialize(G4int);
  G4double GetOutOfWorldTracks(G4int, G4int);
  G4double GetSecondaries(G4int, G4int);
  G4double GetEdepFCAL(G4String);

};

#endif
