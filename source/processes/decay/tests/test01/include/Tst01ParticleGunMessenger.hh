//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: Tst01ParticleGunMessenger.hh,v 1.2 2001-07-11 10:02:28 gunter Exp $
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








