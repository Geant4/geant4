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
// $Id: Tst23PhysicsListMessenger.hh,v 1.1 2001-12-14 14:53:24 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst23PhysicsListMessenger_h
#define Tst23PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst23PhysicsList;
class G4UIcommand; 
class G4UIdirectory;

class Tst23PhysicsListMessenger: public G4UImessenger
{
  public:
    Tst23PhysicsListMessenger(Tst23PhysicsList* myPL);
    virtual ~Tst23PhysicsListMessenger();
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String  GetCurrentValue(G4UIcommand * command);

  private:
    Tst23PhysicsList*           myPhysList;
    G4UIcommand *                cmdSetCut;
    G4UIdirectory*               myphysDir;
};

#endif

