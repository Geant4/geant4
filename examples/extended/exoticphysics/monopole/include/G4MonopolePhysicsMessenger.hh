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
// $Id: G4MonopolePhysicsMessenger.hh,v 1.2 2010-11-29 15:14:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  12.07.10  S.Burdin (changed the magnetic and electric charge variables from integer to double)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4MonopolePhysicsMessenger_h
#define G4MonopolePhysicsMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4MonopolePhysics;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MonopolePhysicsMessenger: public G4UImessenger
{
public:

  G4MonopolePhysicsMessenger(G4MonopolePhysics*);
  ~G4MonopolePhysicsMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:

  G4MonopolePhysics*    phys;
    
  G4UIdirectory*             mPhysicsDir;    
  G4UIcommand*               mPhysicsCmd;
  G4UIcmdWithADouble*      mCmd;
  G4UIcmdWithADouble*      zCmd;
  G4UIcmdWithADoubleAndUnit* massCmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
