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
// $Id: EventActionMessenger.hh,v 1.1 2003-07-31 01:16:04 dwright Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#ifndef EventActionMessenger_h
#define EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class EventAction;
class G4UIcmdWithAnInteger;


class EventActionMessenger: public G4UImessenger
{
  public:
    EventActionMessenger(EventAction*);
    ~EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    EventAction*   eventAction;   
    G4UIcmdWithAnInteger* PrintCmd;    
};

#endif
