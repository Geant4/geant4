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
//
// $Id: G3EleTable.hh,v 1.7 2001-07-16 15:38:18 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

extern G3EleTable G3Ele;
#endif
