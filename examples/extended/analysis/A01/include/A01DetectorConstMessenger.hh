// $Id: A01DetectorConstMessenger.hh,v 1.1 2002-11-13 07:17:38 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#ifndef A01DetectorConstMessenger_h
#define A01DetectorConstMessenger_h 1

class A01DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"

class A01DetectorConstMessenger: public G4UImessenger
{
  public:
    A01DetectorConstMessenger(A01DetectorConstruction* mpga);
    ~A01DetectorConstMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    A01DetectorConstruction * target;

  private: //commands
    G4UIdirectory *             mydetDirectory;
    G4UIcmdWithADoubleAndUnit*  armCmd;

};

#endif


