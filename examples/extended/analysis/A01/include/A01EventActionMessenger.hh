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
// $Id: A01EventActionMessenger.hh,v 1.3 2002-12-13 11:34:28 gunter Exp $
// --------------------------------------------------------------
//
#ifndef A01EventActionMessenger_h
#define A01EventActionMessenger_h 1

class A01EventAction;
class G4UIcmdWithAnInteger;

#include "G4UImessenger.hh"
#include "globals.hh"

class A01EventActionMessenger: public G4UImessenger
{
  public:
    A01EventActionMessenger(A01EventAction* mpga);
    ~A01EventActionMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    A01EventAction* target;

  private: //commands
    G4UIcmdWithAnInteger*  verboseCmd;

};

#endif


