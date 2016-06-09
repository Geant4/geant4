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
#ifndef exrdmMaterialMessenger_h
#define exrdmMaterialMessenger_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4UImessenger.hh"

#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

class exrdmMaterial;
////////////////////////////////////////////////////////////////////////////////
//
class exrdmMaterialMessenger: public G4UImessenger
{
public:
  exrdmMaterialMessenger(exrdmMaterial* );
  ~exrdmMaterialMessenger();

  void SetNewValue (G4UIcommand*, G4String);

private:

  exrdmMaterial                *materialsManager;

  G4UIdirectory             *MaterialDir;
  G4UIcmdWithoutParameter   *ListCmd;
  G4UIcmdWithAnInteger      *DeleteIntCmd;
  G4UIcmdWithAString        *DeleteNameCmd;
  G4UIcommand               *AddCmd;
};
////////////////////////////////////////////////////////////////////////////////
#endif
