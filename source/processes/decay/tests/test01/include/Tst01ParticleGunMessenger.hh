// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst01ParticleGunMessenger.hh,v 1.1 2001-02-08 08:41:44 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef Tst01ParticleGunMessenger_h
#define Tst01ParticleGunMessenger_h 1

class Tst01ParticleGun;

class G4ParticleTable;
class G4UIcommand;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"


class Tst01ParticleGunMessenger: public  G4UImessenger
{
  public:
    Tst01ParticleGunMessenger(Tst01ParticleGun * fPtclGun);
    ~Tst01ParticleGunMessenger();
    
  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    Tst01ParticleGun*   fParticleGun;
    G4ParticleTable*    particleTable;

  private: //commands
    G4UIcmdWithADoubleAndUnit * properTimeCmd;
  
};

#endif








