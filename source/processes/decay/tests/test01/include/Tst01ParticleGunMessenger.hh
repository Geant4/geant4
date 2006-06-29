//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: Tst01ParticleGunMessenger.hh,v 1.3 2006-06-29 19:31:50 gunter Exp $
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








