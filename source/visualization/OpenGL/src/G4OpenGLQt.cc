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
// $Id:$
// GEANT4 tag $Name: not supported by cvs2svn $
//
// John Allison  27th October 2012
// Base class for OpenGLImmediate/StoredQt graphics system factories.

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#include "G4OpenGLQt.hh"

#include "G4UIQt.hh"
#include "G4UImanager.hh"
#include "G4UIbatch.hh"

G4OpenGLQt::G4OpenGLQt (const G4String& name,
                        const G4String& nickname,
                        const G4String& description,
                        Functionality f):
G4VGraphicsSystem (name,
                   nickname,
                   description,
                   f)
{}

G4bool G4OpenGLQt::IsUISessionCompatible () const
{
  G4bool isCompatible = false;
  G4UImanager* ui = G4UImanager::GetUIpointer();
  G4UIsession* session = ui->GetSession();
  // In case it's a G4UIbatch, find original session by recursive search until
  // the session is no longer a G4UIbatch, in which case it will be the
  // original session, if any.
  while (G4UIbatch* batch = dynamic_cast<G4UIbatch*>(session)) {
    session = batch->GetPreviousSession();
  }
  if (!session) {
    // The user has not instantiated a session - must be batch.
    // It's OK to have a Qt window in batch - you can open a viewer, create a
    // scene, set view parameters and /vis/ogl/printEPS.
    isCompatible = true;
  } else {
    // The user has instantiated a session...
    if (dynamic_cast<G4UIQt*>(session)) {
      // ...and it's a G4UIQt session, which is OK.
      isCompatible = true;
    } else {
      // Not OK - go and find the fallback graphics system.
      isCompatible = false;
    }
  }
  return isCompatible;
}

#endif
