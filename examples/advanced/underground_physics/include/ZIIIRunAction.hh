// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ZIIIRunAction.hh,v 1.1 2001-06-26 11:23:21 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef ZIIIRunAction_h
#define ZIIIRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ZIIIRunActionMessenger;
class G4Run;

class ZIIIRunAction : public G4UserRunAction
{
  public:
    ZIIIRunAction();
   ~ZIIIRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

    void SetFilename(G4String val) {fileName = val;};

  private:
  
  G4String fileName;  // log file name for the run
  ZIIIRunActionMessenger* runMessenger;



};

#endif

