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
// $Id: G4UIcsh.hh,v 1.3 2001-07-11 10:01:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4UIcsh_h
#define G4UIcsh_h

#include "G4VUIshell.hh"

//
//   Description:
//   This class gives csh-like shell. (same as dumb terminal)
//

class G4UIcsh : public G4VUIshell {
protected:

public:
  G4UIcsh(const G4String& prompt="%s> ");
  ~G4UIcsh();
  
  virtual G4String GetCommandLine(const char* msg=0);

};

#endif

