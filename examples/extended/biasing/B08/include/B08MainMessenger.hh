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
// $Id: B08MainMessenger.hh,v 1.1 2002/06/04 11:14:51 dressel Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

#ifndef B08MainMessenger_hh
#define B08MainMessenger_hh B08MainMessenger_hh

#include "G4UImessenger.hh"
#include "globals.hh"
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcommand;


class B08MainMessenger: public G4UImessenger{
public: 
  B08MainMessenger();
  ~B08MainMessenger(){}
  void SetNewValue(G4UIcommand* pCmd,G4String szValue);
  G4String GetResultsFileName() const {return fResultsFileName;}
  G4int GetNumberOfEvents() const {return fNumberOfEvents;}
  G4double GetUpperLimit() const {return fUpper;}
  G4double GetLowerLimit() const {return fLower;}

private:
  G4String fResultsFileName;
  G4int fNumberOfEvents;
  G4double fUpper;
  G4double fLower;

  G4UIcmdWithAString *fCmdFn;
  G4UIcmdWithAnInteger *fCmdEv;
  G4UIcmdWithADouble *fCmdUpper;
  G4UIcmdWithADouble *fCmdLower;
  

};


#endif
