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
// $Id: G4HepRepViewer.cc,v 1.1 2001-08-24 23:06:42 perl Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4HepRepViewer.hh"

#include "G4ios.hh"
#include "g4std/strstream"

#include "G4VSceneHandler.hh"

#include "G4HepRepSceneHandler.hh"

//HepRep
#include "JHepRep.hh"
#include "JHepRepFactory.hh"

G4HepRepViewer::G4HepRepViewer
(G4VSceneHandler& sceneHandler, const G4String& name):
  G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name) {

  factory = ((G4HepRepSceneHandler*)(&sceneHandler))->GetHepRepFactory();
  heprep = ((G4HepRepSceneHandler*)(&sceneHandler))->GetHepRep();
}

G4HepRepViewer::~G4HepRepViewer() {}

void G4HepRepViewer::SetView() {
  G4cout << "G4HepRepViewer::SetView() called." << G4endl;
}

void G4HepRepViewer::ClearView() {
  G4cout << "G4HepRepViewer::ClearView() called." << G4endl;
}

void G4HepRepViewer::DrawView() {
  G4cout << "G4HepRepViewer::DrawView() called." << G4endl;
  NeedKernelVisit ();  // Always need to visit G4 kernel.
  ProcessView ();
}

void G4HepRepViewer::ShowView () {
//    OCameraViewNode (fGoCamera,fSceneHandler.GetRootNode());  
#ifdef DEBUG
    G4cout << "G4HepRepViewer::ShowView" << G4endl;
#endif
    G4VViewer::ShowView();
    
    if (factory->SaveAsXML(heprep, "HepRepTest.xml")) {
        // FIXME: no need to exit JVM file here
        factory->Error("Could not write XML file");
    }
}
