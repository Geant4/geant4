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
// $Id: HadrontherapyPrimaryGeneratorMessenger.hh,v 1.2 2005-04-28 20:39:33 mpiergen Exp $
//
//
// Code developed by: M. Piergentili

#ifndef HadrontherapyPrimaryGeneratorMessenger_h
#define HadrontherapyPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class HadrontherapyPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

//**********************************************************************

class HadrontherapyPrimaryGeneratorMessenger: public G4UImessenger
{
  public:

    HadrontherapyPrimaryGeneratorMessenger(HadrontherapyPrimaryGeneratorAction*);
   ~HadrontherapyPrimaryGeneratorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);

  private:
    HadrontherapyPrimaryGeneratorAction* HadrontherapyAction;
    G4UIdirectory*              beamDir;
    G4UIcmdWithADoubleAndUnit*  EnergyCmd;
    G4UIcmdWithADoubleAndUnit*  SourceTypeCmd;
    G4UIcmdWithADoubleAndUnit*  BeamCmd;
    G4UIcmdWithADoubleAndUnit*  SizeYCmd;
    G4UIcmdWithADoubleAndUnit*  SizeZCmd;
};

//**********************************************************************

#endif

