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
//
// John Allison  16th July 2015
// Base class for OpenGLImmediate/StoredXm graphics system factories.

#include "G4OpenGLXm.hh"

#include "G4UImanager.hh"
#include "G4UIbatch.hh"

G4OpenGLXm::G4OpenGLXm (const G4String& name,
                        const G4String& nickname,
                        const G4String& description,
                        Functionality f):
G4VGraphicsSystem (name,
                   nickname,
                   description,
                   f)
{}

G4bool G4OpenGLXm::IsUISessionCompatible () const
{
  // Xm windows are not appropriate in a batch session.
  G4UIsession* baseSession = G4UImanager::GetUIpointer()->GetBaseSession();
  if (baseSession == nullptr  // Pure batch session
      || dynamic_cast<G4UIbatch*>(baseSession) != nullptr) return false;
  return true;
}
