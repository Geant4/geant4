// $Id: A01MagneticFieldMessenger.hh,v 1.1 2002-11-13 07:18:42 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#ifndef A01MagneticFieldMessenger_h
#define A01MagneticFieldMessenger_h 1

class A01MagneticField;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"

class A01MagneticFieldMessenger: public G4UImessenger
{
  public:
    A01MagneticFieldMessenger(A01MagneticField* mpga);
    ~A01MagneticFieldMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    A01MagneticField * target;

    G4UIcmdWithADoubleAndUnit*  fieldCmd;

};

#endif


