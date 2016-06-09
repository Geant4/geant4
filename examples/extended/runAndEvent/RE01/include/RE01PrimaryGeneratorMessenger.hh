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
// $Id: RE01PrimaryGeneratorMessenger.hh,v 1.1 2004/11/26 07:37:40 asaim Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//


#ifndef RE01PrimaryGeneratorMessenger_h
#define RE01PrimaryGeneratorMessenger_h 1

class RE01PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;

#include "G4UImessenger.hh"
#include "globals.hh"

class RE01PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    RE01PrimaryGeneratorMessenger(RE01PrimaryGeneratorAction* mpga);
    ~RE01PrimaryGeneratorMessenger();
    
  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    RE01PrimaryGeneratorAction * myAction;
    
  private: //commands
    G4UIdirectory *             mydetDirectory;
    G4UIcmdWithAString *        genCmd;
    
};

#endif


