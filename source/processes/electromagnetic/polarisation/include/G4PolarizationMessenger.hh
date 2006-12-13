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
// -------------------------------------------------------------------
//
// $Id: G4PolarizationMessenger.hh,v 1.2 2006-12-13 15:44:33 gunter Exp $
// tag $Name: not supported by cvs2svn $
//
// GEANT4 Class header file
//
//
// File name:     G4PolarizationMessenger
//
// Author:        Andreas Schaelicke
//
// Creation date: 01.05.2005
//
// Modifications:
//
// Class Description:
//
// Provides access to general polarization information and to 
// polarization for logical volumes through macro files.
//

#ifndef G4PolarizationMessenger_h
#define G4PolarizationMessenger_h 1

class G4PolarizationManager;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4UIcommand;

#include "G4UImessenger.hh"
#include "globals.hh"

class G4PolarizationMessenger: public G4UImessenger
{
  public:
    G4PolarizationMessenger(G4PolarizationManager* runMgr);
    ~G4PolarizationMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);
  private:
    G4PolarizationManager * polarizationManager;
    
  private: //commands
    G4UIdirectory *             polarizationDirectory;

    G4UIdirectory *             managerDirectory;
    G4UIcmdWithAnInteger *      verboseCmd;
    G4UIcmdWithABool *          optActivateCmd;
    
    G4UIdirectory *             volumeDirectory;
    G4UIcmdWithoutParameter *   printVolumeListCmd;
    G4UIcommand *               setPolarizationCmd;

    G4UIdirectory *             testDirectory;
    G4UIcmdWithoutParameter *   testPolarizationTransformationCmd;
    G4UIcmdWithoutParameter *   testInteractionFrameCmd;
};

#endif


