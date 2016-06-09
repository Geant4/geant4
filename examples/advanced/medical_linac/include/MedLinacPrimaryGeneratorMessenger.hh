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
// $Id: MedLinacPrimaryGeneratorMessenger.hh,v 1.3 2004/06/18 09:17:40 gunter Exp $
//
//
// Code developed by: M. Piergentili

#ifndef MedLinacPrimaryGeneratorMessenger_h
#define MedLinacPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class MedLinacPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

//**********************************************************************

class MedLinacPrimaryGeneratorMessenger: public G4UImessenger
{
  public:

    MedLinacPrimaryGeneratorMessenger(MedLinacPrimaryGeneratorAction*);
   ~MedLinacPrimaryGeneratorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);

  private:
    MedLinacPrimaryGeneratorAction* MedLinacAction;
    G4UIcmdWithADoubleAndUnit*  EnergyCmd;
    G4UIcmdWithADoubleAndUnit*  SourceTypeCmd;
};

//**********************************************************************

#endif

