//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4NeutronKillerMessenger.hh 68048 2013-03-13 14:34:07Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4NeutronKillerMessenger
//
// Description: Messenger class
//
// Author:      V.Ivanchenko 26/09/00 for HARP software
//
//----------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4NeutronKillerMessenger_h
#define G4NeutronKillerMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4NeutronKiller;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4NeutronKillerMessenger: public G4UImessenger
{
public:

  G4NeutronKillerMessenger(G4NeutronKiller*);
  virtual ~G4NeutronKillerMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:

  // hide assignment operator as private
  G4NeutronKillerMessenger(const G4NeutronKillerMessenger&);
  G4NeutronKillerMessenger& operator = (const G4NeutronKillerMessenger &right);

  G4NeutronKiller*   killer;
    
  G4UIdirectory* dir;
  G4UIcmdWithADoubleAndUnit* eCmd;
  G4UIcmdWithADoubleAndUnit* tCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
