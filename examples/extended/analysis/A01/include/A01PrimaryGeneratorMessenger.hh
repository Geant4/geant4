// $Id: A01PrimaryGeneratorMessenger.hh,v 1.1 2002-11-13 07:18:56 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#ifndef A01PrimaryGeneratorMessenger_h
#define A01PrimaryGeneratorMessenger_h 1

class A01PrimaryGeneratorAction;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;

#include "G4UImessenger.hh"
#include "globals.hh"

class A01PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    A01PrimaryGeneratorMessenger(A01PrimaryGeneratorAction* mpga);
    ~A01PrimaryGeneratorMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    A01PrimaryGeneratorAction * target;

  private: //commands
    G4UIcmdWithADoubleAndUnit*  momentumCmd;
    G4UIcmdWithADoubleAndUnit*  sigmaMomCmd;
    G4UIcmdWithADoubleAndUnit*  sigmaAngCmd;
    G4UIcmdWithABool*           randomCmd;

};

#endif


