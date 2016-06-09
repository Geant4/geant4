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
//    ************************************************
//    *                                              *
//    *      RadmonGunMessenger.hh      *
//    *                                              *
//    ************************************************
//
// $Id: RadmonGunMessenger.hh,v 1.2.2.2 2009/08/11 14:20:35 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-04 $
//
//Code developed by:  S.Guatelli, guatelli@ge.infn.it
//
 
#ifndef RadmonGunMessenger_h
#define RadmonGunMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class RadmonPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;

class RadmonGunMessenger: public G4UImessenger
{
public:
  RadmonGunMessenger(RadmonPrimaryGeneratorAction* );
  ~RadmonGunMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  RadmonPrimaryGeneratorAction* primary;
  G4UIdirectory*                gunDir; 
  G4UIcmdWithAString*           dataCmd;
};
#endif

