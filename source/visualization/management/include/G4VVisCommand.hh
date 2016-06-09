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
// $Id: G4VVisCommand.hh,v 1.17 2005/03/09 23:48:15 allison Exp $
// GEANT4 tag $Name: geant4-07-01 $

// Base class for visualization commands - John Allison  9th August 1998
// It is really a messenger - we have one command per messenger.

#ifndef G4VVISCOMMAND_HH
#define G4VVISCOMMAND_HH

#include "G4VisManager.hh"
#include "G4UImessenger.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4UIcommand;
class G4UIcmdWithAString;

class G4VVisCommand: public G4UImessenger {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VVisCommand ();
  virtual ~G4VVisCommand ();
  static void SetVisManager (G4VisManager*);

protected:

  // Conversion routines augmenting those in G4UIcommand.
  static G4String ConvertToString(G4double x, G4double y,
				  const char * unitName);
  static void  ConvertToDoublePair(const G4String& paramString,
				   G4double& xval,
				   G4double& yval);

  // Other utilities.
  void UpdateVisManagerScene (const G4String& sceneName = "");

  // Data members.
  static G4VisManager* fpVisManager;
};

#include "G4VVisCommand.icc"

#endif
