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
///////////////////////////////////////////////////////////////////////////////
// File: CCalSensAssign.hh
// Description: CCalSenAssign creates and assigns the sensitive detetctors 
//              from the map of logical volumes which are potentially sensitive
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalSensAssign_h
#define CCalSensAssign_h

#include "CCalVOrganization.hh"

#include "g4std/map"
#include "globals.hh"
#include "G4VSensitiveDetector.hh"

class CCalSensAssign {

public:
  ~CCalSensAssign(){};
  static CCalSensAssign* getInstance();
  bool assign();
  bool stackingAction();
  bool addCaloSD(G4String name, CCalVOrganization* numberingScheme);

private:
  CCalSensAssign();

  static CCalSensAssign* theInstance;
  G4std::map<G4String,G4VSensitiveDetector*> sens_;

};
#endif
