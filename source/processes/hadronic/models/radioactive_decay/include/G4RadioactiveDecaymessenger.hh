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
#ifndef G4RadioactiveDecaymessenger_h
#define G4RadioactiveDecaymessenger_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              RadioactiveDecaymessenger.hh
//
// Version:             0.b.4
// Date:                14/04/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// 13 April 2000, F Lei, DERA UK
// 0.b.4 release. No change to this file     
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
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

#include "G4RadioactiveDecay.hh"
#include "G4UIcmdWithNucleusLimits.hh"

class G4RadioactiveDecay;
////////////////////////////////////////////////////////////////////////////////
//
class G4RadioactiveDecaymessenger: public G4UImessenger
{
  // class description
  // The G4RadioactiveDecaymessenger is instantiated by the G4RadioactiveDecay
  // process and introduces into the UI menu commands to control the running
  // of G4RadioactiveDecay

public: //with description
    G4RadioactiveDecaymessenger (G4RadioactiveDecay*
      theRadioactiveDecayContainer);
  //    Constructor introduces commands into the UI menu to control
  //    G4RadioactiveDecay.  theRadioactiveDecayContainer1 is used to identify
  //    to this class (when instatiated) the associated G4RadioactiveDecay
  //    process whose parameters are going to be changed as a result
  //    of the UI commands.
  ~G4RadioactiveDecaymessenger ();
  //    Destructor deletes G4UIdirectory and G4UIcommand objects.
  //
  void SetNewValue (G4UIcommand *command, G4String newValues);
  //    Identifies the command which has been invoked by the user, extracts the
  //    parameters associated with that command (held in newValues, and uses
  //    these values with the appropriate member function of G4RadioactiveDecay.
  //
 private:
  G4RadioactiveDecay             *theRadioactiveDecayContainer;
  
  G4UIdirectory                  *grdmDirectory;
  G4UIcmdWithNucleusLimits       *nucleuslimitsCmd;
  G4UIcmdWithAString             *sourcetimeprofileCmd;
  G4UIcmdWithAString             *decaybiasprofileCmd;
  G4UIcmdWithABool               *analoguemcCmd;
  G4UIcmdWithABool               *fbetaCmd;
  G4UIcmdWithABool               *brbiasCmd;
  G4UIcmdWithAnInteger           *splitnucleiCmd;
  G4UIcmdWithAnInteger           *verboseCmd;
  G4UIcmdWithAString             *avolumeCmd;
  G4UIcmdWithAString             *deavolumeCmd;
  G4UIcmdWithoutParameter        *allvolumesCmd;
  G4UIcmdWithoutParameter        *deallvolumesCmd;
  G4UIcmdWithABool               *icmCmd;
  G4UIcmdWithABool               *armCmd;
  G4UIcmdWithADoubleAndUnit      *hlthCmd;

  G4UIcommand					 *userDecayDataCmd;
  G4UIcommand					 *userEvaporationDataCmd;

  G4UIcmdWith3Vector             *colldirCmd;
  G4UIcmdWithADoubleAndUnit      *collangleCmd;

};

#endif

