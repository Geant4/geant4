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
// $Id: ExN05EventActionMessenger.hh,v 1.4 2002-01-09 17:24:18 ranjard Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef ExN05EventActionMessenger_h
#define ExN05EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
class G4UIdirectory;
class G4UIcmdWithABool;

class ExN05EventAction;

class ExN05EventActionMessenger: public G4UImessenger
{
  public:
    ExN05EventActionMessenger(ExN05EventAction* EA);
    void SetNewValue(G4UIcommand* command, G4String newValues);
  private:
    ExN05EventAction* EventAction;
    G4UIdirectory*    eventDirectory;
    G4UIcmdWithABool* drawEventCmd;
};

#endif

