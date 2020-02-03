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
///////////////////////////////////////////////////////////////////////////////
// File: CCalSensAssign.hh
// Description: CCalSenAssign creates and assigns the sensitive detetctors 
//              from the map of logical volumes which are potentially sensitive
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalSensAssign_h
#define CCalSensAssign_h 1

#include "CCalVOrganization.hh"

#include <map>
#include "globals.hh"
#include "G4VSensitiveDetector.hh"

class CCalSensAssign
{
public:
  ~CCalSensAssign(){}
  static CCalSensAssign* getInstance();
  G4bool assign();
  G4bool stackingAction();
  G4bool addCaloSD(G4String name, CCalVOrganization* numberingScheme);

private:
  CCalSensAssign();

  static CCalSensAssign* theInstance;
  std::map<G4String,G4VSensitiveDetector*> sens_;

};
#endif
