//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FCALSteppingAction.hh,v 1.4 2002-12-17 15:53:22 pmendez Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
