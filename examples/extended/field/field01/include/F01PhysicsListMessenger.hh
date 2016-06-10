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
/// \file field/field01/include/F01PhysicsListMessenger.hh
/// \brief Definition of the F01PhysicsListMessenger class
//
//
// $Id: F01PhysicsListMessenger.hh 76248 2013-11-08 11:19:52Z gcosmo $
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F01PhysicsListMessenger_h
#define F01PhysicsListMessenger_h 1

#include "G4UImessenger.hh"

class F01PhysicsList;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class F01PhysicsListMessenger: public G4UImessenger
{
  public:
    F01PhysicsListMessenger(F01PhysicsList*);
    virtual ~F01PhysicsListMessenger();
 
    virtual void SetNewValue(G4UIcommand*, G4String);

  private:
    F01PhysicsList*          fF01List;

    G4UIcmdWithADoubleAndUnit* fSetMaxStepCmd;

    G4UIcmdWithADoubleAndUnit* fCutGCmd;
    G4UIcmdWithADoubleAndUnit* fCutECmd;
};

#endif
