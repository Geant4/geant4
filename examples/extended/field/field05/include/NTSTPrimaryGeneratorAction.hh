// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: NTSTPrimaryGeneratorAction instantiates each of the allowed
// generators and invokes the one specified by the messenger

// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef NTSTPrimaryGeneratorAction_hh
#define NTSTPrimaryGeneratorAction_hh 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "globals.hh"

class G4Event;
class NTSTPrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class NTSTPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  NTSTPrimaryGeneratorAction();
  virtual ~NTSTPrimaryGeneratorAction();
  void GeneratePrimaries(G4Event* event);
  void Print(const G4Event*) const;

private:
  NTSTPrimaryGeneratorMessenger* messenger;
};

#endif




