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
/// \file persistency/gdml/G03/include/G03DetectorMessenger.hh
/// \brief Definition of the G03DetectorMessenger class
//
//
//
// Class G03DetectorMessenger
//
// Utility messenger for defining run-time commands relative to the example.
//
// ----------------------------------------------------------------------------

#ifndef G03DetectorMessenger_h
#define G03DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G03DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

/// Detector messenger for the GDML extensions example

class G03DetectorMessenger: public G4UImessenger
{

  public:

    G03DetectorMessenger( G03DetectorConstruction* );
   ~G03DetectorMessenger();
    
    virtual void SetNewValue( G4UIcommand*, G4String );

  private:

    G03DetectorConstruction      *fTheDetector;
    G4UIdirectory             *fTheDetectorDir;
    G4UIcmdWithAString        *fTheReadCommand, *fTheWriteCommand;
};

// ----------------------------------------------------------------------------

#endif
