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
// $Id: G4HepRepFileViewer.cc,v 1.4 2002-01-29 01:25:22 perl Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4HepRepFileViewer.hh"

#include "G4ios.hh"
#include "g4std/strstream"

#include "G4VSceneHandler.hh"

#include "G4HepRepFileSceneHandler.hh"

//HepRep
#include "HepRepXMLWriter.hh"

G4HepRepFileViewer::G4HepRepFileViewer
(G4VSceneHandler& sceneHandler, const G4String& name):
  G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name) {
  fileCounter = 0;
  isOpen = false;
  hepRepXMLWriter = ((G4HepRepFileSceneHandler*)(&sceneHandler))->GetHepRepXMLWriter();
}

G4HepRepFileViewer::~G4HepRepFileViewer() {}

void G4HepRepFileViewer::SetView() {
#ifdef G4HEPREPFILEDEBUG
  G4cout << "G4HepRepFileViewer::SetView() called with isOpen=" << isOpen << G4endl;
#endif
  
  if (isOpen) {    
    hepRepXMLWriter->endTypes();
  } else {
    isOpen = true;
    char* newFileSpec;
    newFileSpec = new char [100];
    int length;
    length = sprintf (newFileSpec, "%s%d%s","G4Data",fileCounter,".heprep");
    hepRepXMLWriter->open(newFileSpec);
    fileCounter++;

    hepRepXMLWriter->addAttDef("LVol", "Logical Volume", "Physics","");
    hepRepXMLWriter->addAttDef("Solid", "Solid Name", "Physics","");
    hepRepXMLWriter->addAttDef("EType", "Entity Type", "Physics","");
    hepRepXMLWriter->addAttDef("Material", "Material Name", "Physics","");
    hepRepXMLWriter->addAttDef("Density", "Material Density", "Physics","");
    hepRepXMLWriter->addAttDef("State", "Material State", "Physics","");
    hepRepXMLWriter->addAttDef("Radlen", "Material Radiation Length", "Physics","");
  }
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
  
  hepRepXMLWriter->close();
  isOpen = false;
}
