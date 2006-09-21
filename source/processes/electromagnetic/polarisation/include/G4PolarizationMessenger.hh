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
// -------------------------------------------------------------------
//
// $Id: G4PolarizationMessenger.hh,v 1.1 2006-09-21 21:35:10 vnivanch Exp $
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


