// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN03RunAction.hh,v 1.1 1999-01-07 16:05:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef ExN03RunAction_h
#define ExN03RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

class ExN03RunAction : public G4UserRunAction
{
  public:
    ExN03RunAction();
   ~ExN03RunAction();

  public:
    void BeginOfRunAction(G4Run*);
    void EndOfRunAction(G4Run*);

  private:
    G4int runIDcounter;
};

#endif

