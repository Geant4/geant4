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
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------

#ifndef Tst50PrimaryGeneratorMessenger_h
#define Tst50PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UIcmdWithAString;
class Tst50PrimaryGeneratorAction;

class Tst50PrimaryGeneratorMessenger: public G4UImessenger
{
 public:
  Tst50PrimaryGeneratorMessenger(Tst50PrimaryGeneratorAction*);
  ~Tst50PrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
 private:
  Tst50PrimaryGeneratorAction* tst50Gun; 
  G4UIcmdWithAString*          randomDirectionCmd;
  //command to choose a primary particle random direction
};

#endif

