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
// $Id: G4EmMessenger.hh,v 1.1 2005-11-11 22:56:07 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmMessenger
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 09.11.2005 V.Ivanchenko edit to provide a standard
//
//----------------------------------------------------------------------------
//

#ifndef G4EmMessenger_h
#define G4EmMessenger_h

class G4EmExtraBuilder;

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

class G4EmMessenger: public G4UImessenger
{
public:
  G4EmMessenger(G4EmExtraBuilder* af);
  virtual ~G4EmMessenger();

  void SetNewValue(G4UIcommand* aComm, G4String aS);

private:
  G4EmExtraBuilder*   theB;
  G4UIcmdWithAString* theSynch;
  G4UIcmdWithAString* theGN;
  G4UIdirectory*      aDir1;
  G4UIdirectory*      aDir2;
};

#endif
