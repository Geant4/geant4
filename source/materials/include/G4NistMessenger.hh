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
// $Id: G4NistMessenger.hh,v 1.1 2005/02/22 10:11:09 maire Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// File name:     G4NistMessenger
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 23.12.2004
//
// Modifications:
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

  G4NistMessenger(G4NistManager* );
 ~G4NistMessenger();

  void SetNewValue(G4UIcommand*, G4String);

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
};

#endif

