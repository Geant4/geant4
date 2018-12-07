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
/// \file persistency/gdml/G02/include/G02DetectorMessenger.hh
/// \brief Definition of the G02DetectorMessenger class
//
//
//
// Class G02DetectorMessenger
//
// Utility messenger for defining run-time commands relative to the example.
//
// ----------------------------------------------------------------------------

#ifndef G02DetectorMessenger_h
#define G02DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G02DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

// ----------------------------------------------------------------------------

/// Detector messenger class used in GDML read/write example

class G02DetectorMessenger: public G4UImessenger
{

  public:

    G02DetectorMessenger( G02DetectorConstruction* );
   ~G02DetectorMessenger();
    
    virtual void SetNewValue( G4UIcommand*, G4String );

  private:

    G02DetectorConstruction*      fTheDetector;
    G4UIdirectory*             fTheDetectorDir;
    G4UIcmdWithAString*        fTheReadCommand;
    G4UIcmdWithAString*        fTheWriteCommand;
    G4UIcmdWithAString*        fTheStepCommand;
};

// ----------------------------------------------------------------------------

#endif
