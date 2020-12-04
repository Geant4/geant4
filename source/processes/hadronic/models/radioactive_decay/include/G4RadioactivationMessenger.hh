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
#ifndef G4RadioactivationMessenger_h
#define G4RadioactivationMessenger_h 1

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4RadioactivationMessenger.hh                                     //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   29 August 2017                                                    //
//  Description: messenger class for biased version of G4RadioactiveDecay.    //
//  Based on the code of F. Lei and P.R. Truscott.                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

#include "G4Radioactivation.hh"
#include "G4UIcmdWithNucleusLimits.hh"

class G4Radioactivation;

class G4RadioactivationMessenger: public G4UImessenger
{
  public:
    G4RadioactivationMessenger(G4Radioactivation* theRadioactivationContainer);
    ~G4RadioactivationMessenger();

    void SetNewValue(G4UIcommand* command, G4String newValues);

  private:
    G4Radioactivation* theRadioactivationContainer;
  
    G4UIdirectory* old_grdmDirectory;              // To be removed in G4 11.0
    G4UIdirectory*      rdmDirectory;            
    G4UIcmdWithABool* old_analoguemcCmd;           // To be removed in G4 11.0
    G4UIcmdWithABool*     analoguemcCmd;
    G4UIcmdWithAString* old_sourcetimeprofileCmd;  // To be removed in G4 11.0
    G4UIcmdWithAString*     sourcetimeprofileCmd;
    G4UIcmdWithAString* old_decaybiasprofileCmd;   // To be removed in G4 11.0
    G4UIcmdWithAString*     decaybiasprofileCmd;
    G4UIcmdWithABool* old_brbiasCmd;               // To be removed in G4 11.0
    G4UIcmdWithABool*     brbiasCmd;
    G4UIcmdWithAnInteger* old_splitnucleiCmd;      // To be removed in G4 11.0
    G4UIcmdWithAnInteger*     splitnucleiCmd;
    G4UIcmdWithADoubleAndUnit* old_hlthCmd;        // To be removed in G4 11.0
    G4UIcmdWithADoubleAndUnit*     hlthCmd;

};

#endif

