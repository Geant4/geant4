// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst18RunAction.hh,v 1.5 2000-06-14 17:48:11 flei Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst18RunAction_h
#define Tst18RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst18RunActionMessenger;
class G4Run;

class Tst18RunAction : public G4UserRunAction
{
  public:
    Tst18RunAction();
   ~Tst18RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
  
    void SetFilename(G4String val) {fileName = val;};

  private:
  
  G4String fileName;  // log file name for the run
  Tst18RunActionMessenger* runMessenger;

};

#endif

