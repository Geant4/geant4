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
#ifndef G4RadioactiveDecayBaseMessenger_h
#define G4RadioactiveDecayBaseMessenger_h 1

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4RadioactiveDecayBaseMessenger.hh                                //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   29 August 2017                                                    //
//  Description: messenger class for the non-biased version of                //
//               G4RadioactiveDecay.  Based on the code of F. Lei and         //
//               P.R. Truscott.                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "globals.hh"

#include "G4RadioactiveDecayBase.hh"
#include "G4UIcmdWithNucleusLimits.hh"

class G4RadioactiveDecayBase;

class G4RadioactiveDecayBaseMessenger: public G4UImessenger
{
  public: //with description
    G4RadioactiveDecayBaseMessenger(G4RadioactiveDecayBase* theRadioactiveDecayContainer);
    ~G4RadioactiveDecayBaseMessenger();

    void SetNewValue (G4UIcommand *command, G4String newValues);

  private:
    G4RadioactiveDecayBase* theRadioactiveDecayContainer;
  
    G4UIdirectory* grdmDirectory;
    G4UIcmdWithNucleusLimits* nucleuslimitsCmd;
    G4UIcmdWithABool* fbetaCmd;
    G4UIcmdWithAnInteger* verboseCmd;
    G4UIcmdWithAString* avolumeCmd;
    G4UIcmdWithAString* deavolumeCmd;
    G4UIcmdWithoutParameter* allvolumesCmd;
    G4UIcmdWithoutParameter* deallvolumesCmd;
    G4UIcmdWithABool* icmCmd;
    G4UIcmdWithABool* armCmd;

    G4UIcommand* userDecayDataCmd;
    G4UIcommand* userEvaporationDataCmd;

    G4UIcmdWith3Vector* colldirCmd;
    G4UIcmdWithADoubleAndUnit* collangleCmd;
};

#endif

