// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TstVAPrimaryGeneratorMessenger.hh,v 1.3 2001-02-07 17:30:59 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------

#ifndef TstVAPrimaryGeneratorMessenger_h
#define TstVAPrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class TstVAPrimaryGeneratorAction;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;

class TstVAPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    TstVAPrimaryGeneratorMessenger(TstVAPrimaryGeneratorAction * myPGA);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    TstVAPrimaryGeneratorAction * myPGAction;
    G4UIcmdWithoutParameter * standardGun;
    G4UIcmdWithAString * randomGun;
};

#endif

