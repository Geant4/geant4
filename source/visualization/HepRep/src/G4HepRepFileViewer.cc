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
// $Id: G4HepRepFileViewer.cc 66373 2012-12-18 09:41:34Z gcosmo $

#include "G4HepRepFileViewer.hh"

#include "G4VSceneHandler.hh"

#include "G4HepRepFileSceneHandler.hh"

//HepRep
#include "G4HepRepFileXMLWriter.hh"

G4HepRepFileViewer::G4HepRepFileViewer
(G4VSceneHandler& sceneHandler, const G4String& name):
  G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name) {
  hepRepXMLWriter = ((G4HepRepFileSceneHandler*)(&sceneHandler))->GetHepRepXMLWriter();
  // Make changes to view parameters for HepRep...
  fVP.SetCulling(false);
  fDefaultVP.SetCulling(false);
}

G4HepRepFileViewer::~G4HepRepFileViewer() {
  ShowView ();
}

void G4HepRepFileViewer::SetView() {
#ifdef G4HEPREPFILEDEBUG
  G4cout << "G4HepRepFileViewer::SetView() called.=" << G4endl;
#endif
}

void G4HepRepFileViewer::ClearView() {
#ifdef G4HEPREPFILEDEBUG
  G4cout << "G4HepRepFileViewer::ClearView() called." << G4endl;
#endif
}

void G4HepRepFileViewer::DrawView() {
#ifdef G4HEPREPFILEDEBUG
  G4cout << "G4HepRepFileViewer::DrawView() called." << G4endl;
#endif
  NeedKernelVisit ();  // Always need to visit G4 kernel.
  ProcessView ();
}

void G4HepRepFileViewer::ShowView () {
#ifdef G4HEPREPFILEDEBUG
  G4cout << "G4HepRepFileViewer::ShowView" << G4endl;
#endif
  G4VViewer::ShowView();
  
  if (hepRepXMLWriter->isOpen)
    hepRepXMLWriter->close();
}
