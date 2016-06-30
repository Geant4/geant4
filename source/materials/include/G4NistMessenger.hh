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
// $Id: G4NistMessenger.hh 96794 2016-05-09 10:09:30Z gcosmo $
//
// File name:     G4NistMessenger
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 23.12.2004
//
// Modifications:
// 29.04.07 V.Ivanchenko add recCmd
//
//
// Class Description:
//  This is a messenger class to interface to exchange information
//  between G4NistManager and UI.
//
//  /material/  directory
//
//   Commands :
//    verbose  val        verbose level for material/element builders
//
//
//  /material/nist/  directory
//
//   Commands :
//    printElement   symbol    output element parameters by symbol
//    printZElement    Z       output element parameters by number
//    listMaterials   key      list of materials (key = simple compound hep)
//
//  /material/g4/  directory
//
//   Commands :
//    printElement   name        output element  parameters
//    printMaterial  name        output material parameters

// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4NistMessenger_h
#define G4NistMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4NistManager;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4NistMessenger: public G4UImessenger
{
public:

  explicit G4NistMessenger(G4NistManager* );
  virtual ~G4NistMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String) final;

private:

  G4NistManager*             manager;
  
  G4UIdirectory*             matDir;
  G4UIcmdWithAnInteger*      verCmd;
    
  G4UIdirectory*             nistDir;
  G4UIcmdWithAString*        prtElmCmd;  
  G4UIcmdWithAnInteger*      przElmCmd;
  G4UIcmdWithAString*        lisMatCmd;
    
  G4UIdirectory*             g4Dir;
  G4UIcmdWithAString*        g4ElmCmd;   
  G4UIcmdWithAString*        g4MatCmd;
  G4UIcmdWithAString*        g4DensCmd;    
};

#endif

