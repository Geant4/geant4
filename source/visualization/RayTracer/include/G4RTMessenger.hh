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
// $Id: G4RTMessenger.hh,v 1.4 2001-07-11 10:09:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

// class description:
//
//  This is a concrete class of G4UImessenger. This class defines various
// UI commands which are unique to G4RayTracer.
//


#ifndef G4RTMessenger_HH
#define G4RTMessenger_HH 1

#include "G4UImessenger.hh"
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4RayTracer;
class G4RTSteppingAction;

class G4RTMessenger : public G4UImessenger
{
  public:
    G4RTMessenger(G4RayTracer* p1,G4RTSteppingAction* p2);
    virtual ~G4RTMessenger();
    
    virtual G4String GetCurrentValue(G4UIcommand * command);
    virtual void SetNewValue(G4UIcommand * command,G4String newValue);

  private:
    G4RayTracer* theTracer;
    G4RTSteppingAction* theSteppingAction;

    G4UIdirectory* rayDirectory;
    G4UIcmdWithAnInteger* columnCmd;
    G4UIcmdWithAnInteger* rowCmd;
    G4UIcmdWith3VectorAndUnit* targetCmd;
    G4UIcmdWith3VectorAndUnit* eyePosCmd;
    G4UIcmdWith3Vector* lightCmd;
    G4UIcmdWithADoubleAndUnit* spanXCmd;
    G4UIcmdWithADoubleAndUnit* headCmd;
    G4UIcmdWithADoubleAndUnit* attCmd;
    G4UIcmdWithABool* distCmd;
    G4UIcmdWithABool* transCmd;
    G4UIcmdWithAString* fileCmd;
};

#endif



