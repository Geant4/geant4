// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: NTSTRunAction.hh,v 1.1 2003-11-07 21:30:28 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef NTSTRunAction_h
#define NTSTRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

class NTSTRunAction : public G4UserRunAction
{
  public:
    NTSTRunAction();
   ~NTSTRunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

  private:
    G4int runIDcounter;
};

#endif

