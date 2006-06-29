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
// $Id: ExN07DetectorMessenger.hh,v 1.6 2006-06-29 17:54:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef ExN07DetectorMessenger_h
#define ExN07DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class ExN07DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;

class ExN07DetectorMessenger: public G4UImessenger
{
  public:
    ExN07DetectorMessenger(ExN07DetectorConstruction* );
    virtual ~ExN07DetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    virtual G4String GetCurrentValue(G4UIcommand * command);
    
  private:
    void UpdateMaterialList();

    ExN07DetectorConstruction* ExN07Detector;
    
    G4UIdirectory*             N07Dir;
    G4UIcmdWithAString*        AbsMaterCmd;
    G4UIcmdWithAString*        GapMaterCmd;
    G4UIcmdWithAnInteger*      numLayerCmd;
    G4UIcmdWithABool*          SerialCmd;
    G4UIcmdWithAnInteger*      verboseCmd;
    G4UIcmdWithABool*          AddMatCmd;
};


#endif

