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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmQMessenger
//
// Author: 16-Oct-2012 A. Ribon
//         Copied from G4EmMessenger and renamed
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4EmQMessenger_h
#define G4EmQMessenger_h 1

class G4EmQExtraPhysics;

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

class G4EmQMessenger: public G4UImessenger
{
public:
  G4EmQMessenger(G4EmQExtraPhysics* af);
  virtual ~G4EmQMessenger();

  void SetNewValue(G4UIcommand* aComm, G4String aS);

private:
  G4EmQExtraPhysics*  theB;
  G4UIcmdWithAString* theSynch;
  G4UIcmdWithAString* theGN;
  G4UIcmdWithAString* theMUN;
  G4UIdirectory*      aDir1;
  G4UIdirectory*      aDir2;
};

#endif
