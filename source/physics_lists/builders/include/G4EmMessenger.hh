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
// $Id: G4EmMessenger.hh,v 1.1 2006-10-31 11:35:01 gunter Exp $
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
// 19.06.2006 V.Ivanchenko add mu-nuclear process
//
//----------------------------------------------------------------------------
//

#ifndef G4EmMessenger_h
#define G4EmMessenger_h 1

class G4EmExtraPhysics;

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

class G4EmMessenger: public G4UImessenger
{
public:
  G4EmMessenger(G4EmExtraPhysics* af);
  virtual ~G4EmMessenger();

  void SetNewValue(G4UIcommand* aComm, G4String aS);

private:
  G4EmExtraPhysics*   theB;
  G4UIcmdWithAString* theSynch;
  G4UIcmdWithAString* theGN;
  G4UIcmdWithAString* theMUN;
  G4UIdirectory*      aDir1;
  G4UIdirectory*      aDir2;
};

#endif
