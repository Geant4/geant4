//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FCALSteppingAction_h
#define FCALSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FCALSteppingAction : public G4UserSteppingAction
{
  public:
    FCALSteppingAction();
    virtual ~FCALSteppingAction();

    virtual void UserSteppingAction(const G4Step*);

  private:
  G4int EventNo;
  G4int IDnow, IDold, IDout;
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

