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
#ifndef G4UIcmdWithNucleusLimits_h
#define G4UIcmdWithNucleusLimits_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4UIcmdWithNucleusAndUnit.hh
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
// Class Description
//
// The G4UIcmdWithNucleusLimits permits input of the range of nuclides
// that will be treated in the G4RadioactiveDecay Class.  Input is expected in
// the form:
//
//			AMin  AMax  ZMin  ZMax
//
// where AMin, AMax, ZMin, ZMax are respectively the upper and lower limits
// in A (atomic mass), and the upper and lower limits in Z (atomic number)
// which define the range of nuclides treated by the G4RadioactiveDecay.
//
// Class Description - end:
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// 
// 13 April 2000, F Lei, DERA UK
// 0.b.4 release. No change to this file        
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "G4UIcommand.hh"
#include "globals.hh"

#include "G4NucleusLimits.hh"
////////////////////////////////////////////////////////////////////////////////
//
class G4UIcmdWithNucleusLimits : public G4UIcommand
{
public: // With description
    G4UIcmdWithNucleusLimits
    (const char * theCommandPath,G4UImessenger * theMessenger);
  //    Constructor identifying the command path in the User Interface and the
  //    associated G4UImessenger which will use this G4UIcommand object.
  //
  ~G4UIcmdWithNucleusLimits();
  //  Destructor
  //
    G4NucleusLimits GetNewNucleusLimitsValue(G4String paramString);
  //    Extracts the values aMin, aMax, zMin, zMax from paramString.
  //    Values returned have aMin, aMax, zMin, and zMax (within the
  //    G4NucleusLimits variable).
    G4String ConvertToString(G4NucleusLimits nuclimit);
  //    Converts the G4NucleusLimits defined by vec into a G4String.
    void SetParameterName(const char * theNameAMin,const char * theNameAMax,
                          const char * theNameZMin, const char * theNameZMax,
                          G4bool omittable, G4bool currentAsDefault=true);
  //    Identifies the parameter names associated with each of the G4doubles
  //    used to define the G4NucleusLimits variable.
    void SetDefaultValue(G4NucleusLimits defVal);
  //    Sets the default G4NucleusLimits if the command is invoked without any
  //    parameters.
};
////////////////////////////////////////////////////////////////////////////////
#endif

