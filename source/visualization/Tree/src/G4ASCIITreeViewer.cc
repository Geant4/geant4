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
// $Id: G4ASCIITreeViewer.cc,v 1.5 2002-12-11 16:08:50 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4ASCIITreeViewer.hh"

#include "G4ios.hh"
#include "g4std/strstream"

G4ASCIITreeViewer::G4ASCIITreeViewer
(G4VSceneHandler& sceneHandler, const G4String& name):
  G4VTreeViewer(sceneHandler, name) {
  // Make changes to view parameters for ASCIITree...
  fVP.SetCulling(false);
  fDefaultVP.SetCulling(false);
}

G4ASCIITreeViewer::~G4ASCIITreeViewer() {}
