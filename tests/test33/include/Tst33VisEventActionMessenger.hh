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
// $Id: Tst33VisEventActionMessenger.hh,v 1.2 2002-10-31 08:32:44 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Tst33VisEventActionMessenger_h
#define Tst33VisEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst33VisEventAction;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Tst33VisEventActionMessenger: public G4UImessenger
{
public:
  explicit Tst33VisEventActionMessenger(Tst33VisEventAction*);
  virtual ~Tst33VisEventActionMessenger();
  
  virtual void SetNewValue(G4UIcommand*, G4String);
  
private:
  Tst33VisEventActionMessenger(const Tst33VisEventActionMessenger &);
  Tst33VisEventActionMessenger &
  operator=(const Tst33VisEventActionMessenger &);
  Tst33VisEventAction*   eventAction;   
  G4UIcmdWithAString* DrawCmd;
  G4UIcmdWithAnInteger* PrintCmd;    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
