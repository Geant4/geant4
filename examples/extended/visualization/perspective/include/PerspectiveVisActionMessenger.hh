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
// $Id: PerspectiveVisActionMessenger.hh,v 1.1 2006-05-25 08:44:52 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef PERSPECTIVEVISACTIONMESSENGER_HH
#define PERSPECTIVEVISACTIONMESSENGER_HH

#include "G4UImessenger.hh"

class PerspectiveVisAction;
class G4UIdirectory;
class G4UIcmdWithAString;

class PerspectiveVisActionMessenger: public G4UImessenger {
public:
  PerspectiveVisActionMessenger(PerspectiveVisAction*);
  ~PerspectiveVisActionMessenger();
  void SetNewValue (G4UIcommand*, G4String);

private:
  PerspectiveVisAction* fPVA;
  G4UIdirectory* fpDirectory;
  G4UIcmdWithAString* fpCommandOS;
  G4UIcmdWithAString* fpCommandScene;
};

#endif

