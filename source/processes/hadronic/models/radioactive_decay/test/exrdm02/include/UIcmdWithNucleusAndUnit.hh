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
#ifndef UIcmdWithNucleusAndUnit_h
#define UIcmdWithNucleusAndUnit_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              UIcmdWithNucleusAndUnit.hh
//
// Version:             0.b.3
// Date:                29/02/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class Description
//
// The G4UIcmdWithNucleusAndUnit permits input of the nucleus definition
// in terms of its (atomic weight, atomic number, excitation energy).
// Input is expected in the form:
//
//			A, Z, E (energy unit)
//
// where A, Z, E are respectively the atomic weight, atomic number and
// excitation energy.  The energy unit can be any of the geant4 defined energy
// units.  The default is "keV"
//
// class description - end
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
//
#include "G4UIcommand.hh"
#include "globals.hh"

#include "Nucleus.hh"
////////////////////////////////////////////////////////////////////////////////
//
class UIcmdWithNucleusAndUnit : public G4UIcommand
{
public:  // with description
  UIcmdWithNucleusAndUnit
  (const char * theCommandPath,G4UImessenger * theMessenger);
  //    Constructor identifying the command path in the User Interface and the
  //    associated G4UImessenger which will use this G4UIcommand object.
  //
  ~UIcmdWithNucleusAndUnit();
  //    Desctructor
  //
  Nucleus GetNewNucleusValue(G4String paramString);
  //    Extracts the values A, Z, E and unit from paramString.
  //
  G4double GetNewUnitValue(G4String paramString);
  //    Returns the value of the unit (paramString) as defined in geant4
  //
  G4String ConvertToString(Nucleus nuc, const char * unit);
  //    Converts the Nuccleus defined by nuc and the associated units of
  //    energy *unit into a G4String.
  void SetParameterName(const char * theNameA,const char * theNameZ,
                        const char * theNameE,
                        G4bool omittable, G4bool currentAsDefault=true);
  //    Identifies the parameter names associated with A, Z, E
  //
  void SetDefaultValue(Nucleus defVal);
  //    Sets the default Nucleus if the command is invoked without any
  //    parameters.
  void SetUnitCandidates(const char * candidateList);
  //    Sets the list of unit candidates
  //
  void SetDefaultUnit(const char * defUnit);
  //    Sets the default unit if the command is invoked without any
  //    parameters.
};
////////////////////////////////////////////////////////////////////////////////
#endif

