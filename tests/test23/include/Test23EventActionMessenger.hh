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
// $Id: Test23EventActionMessenger.hh,v 1.1 2004-03-18 11:02:25 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- Test23TrackingActionMessenger header ----------------
//                 by Mikhail Kossov, December 2003.
//  Test23TrackingActionMessenger class of CHIPS Test of G4QCaptureAtRest process in GEANT4
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------

#ifndef Test23EventActionMessenger_h
#define Test23EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Test23EventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class Test23EventActionMessenger: public G4UImessenger
{
public:
  Test23EventActionMessenger(Test23EventAction*);
  ~Test23EventActionMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:

  Test23EventAction* eventAction;   
  G4UIcmdWithAString* DrawCmd;
  G4UIcmdWithAnInteger* PrintCmd;    
};

#endif
