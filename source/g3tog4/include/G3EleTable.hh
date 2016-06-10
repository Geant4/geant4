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
//
// $Id: G3EleTable.hh 67982 2013-03-13 10:36:03Z gcosmo $
//
// ----------------------
// Class description:
//
// The table of elements.
// In the constructor an array of strings with element name,
// symbol and A are created. 
// The G4Element instances of given Z are created when
// GetEle(G4double Z) is called using the string array.
// For each Z the G4Element is created only once.

// ----------------------

#ifndef G4ELETABLE_HH
#define G4ELETABLE_HH 1

#include "G3toG4Defs.hh"
#include "globals.hh"
#include "G4Element.hh"

class G3EleTable
{

public:  // with description

  G3EleTable();
  virtual ~G3EleTable();
  G4Element* GetEle(G4double Z);

private:

  void LoadUp();
  G4int parse(G4double& Z, char* name, char* sym, G4double& A); 

private:

  char** _EleNames;
  G4Element** _Ele;
  G4int _MaxEle;

};

extern G3G4DLL_API G3EleTable G3Ele;
#endif
