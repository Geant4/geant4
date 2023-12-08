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
//
//---------------------------------------------------------------------------
//
// ClassName:   G4NeutrinoPhysicsMessenger
//
// Author: 2023 V. Ivanchenko created using G4EmMessenger
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4NeutrinoPhysicsMessenger_h
#define G4NeutrinoPhysicsMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"

class G4NeutrinoPhysics;

class G4NeutrinoPhysicsMessenger: public G4UImessenger
{
public:
  explicit G4NeutrinoPhysicsMessenger(G4NeutrinoPhysics* af);
  ~G4NeutrinoPhysicsMessenger() override;

  void SetNewValue(G4UIcommand* aComm, G4String aS) override;

  G4NeutrinoPhysicsMessenger& operator=
  (const G4NeutrinoPhysicsMessenger& right) = delete;
  G4NeutrinoPhysicsMessenger(const G4NeutrinoPhysicsMessenger&) = delete;

private:

  G4NeutrinoPhysics* theB;

  G4UIcmdWithABool* theNu;
  G4UIcmdWithABool* theNuETX;

  G4UIcmdWithADouble* theNuEleCcBF;
  G4UIcmdWithADouble* theNuEleNcBF;
  G4UIcmdWithADouble* theNuNucleusBF;
  G4UIcmdWithADouble* theNuOscDistanceBF;

  G4UIcmdWithAString* theNuDN;
  G4UIcmdWithAString* theNuODN;

  G4UIdirectory* aDir;
};

#endif
