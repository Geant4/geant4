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
//    *      RemSimPrimaryGeneratorMessenger.hh      *
//    *                                              *
//    ************************************************
//
// $Id: RemSimPrimaryGeneratorMessenger.hh,v 1.4 2004/05/22 12:57:05 guatelli Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
//Code developed by:  S.Guatelli, guatelli@ge.infn.it
//
 
#ifndef RemSimPrimaryGeneratorMessenger_h
#define RemSimPrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class RemSimPrimaryGeneratorAction;
class RemSimRunAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class RemSimPrimaryGeneratorMessenger: public G4UImessenger
{
public:
  RemSimPrimaryGeneratorMessenger(RemSimPrimaryGeneratorAction* );
  ~RemSimPrimaryGeneratorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  RemSimPrimaryGeneratorAction* primary;//pointer to detector
  G4UIdirectory*                gunDir; 
  G4UIcmdWithAString*           fluxCmd; //change vehicle 
};
#endif

