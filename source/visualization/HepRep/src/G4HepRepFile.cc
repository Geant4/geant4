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
// $Id: G4HepRepFile.cc,v 1.1 2001-11-06 11:48:07 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Joseph Perl  1 October 2001
// A graphics system to dump geometry hierarchy to the
// HepRep graphics format (HepRep version 1).

//HepRep
#include "HepRepXMLWriter.hh"

#include "G4HepRepFile.hh"
#include "G4HepRepFileSceneHandler.hh"
#include "G4HepRepFileViewer.hh"

static HepRepXMLWriter* hepRepXMLWriter;

G4HepRepFile::G4HepRepFile():
  G4VGraphicsSystem("G4HepRepFile",
		    "HepRepFile",
		    "A template graphics driver",
		    G4VGraphicsSystem::noFunctionality) {
        hepRepXMLWriter = new HepRepXMLWriter();
}

G4HepRepFile::~G4HepRepFile() {}

G4VSceneHandler* G4HepRepFile::CreateSceneHandler(const G4String& name) {
  G4VSceneHandler* pScene = new G4HepRepFileSceneHandler(*this, name);
  G4cout << G4HepRepFileSceneHandler::GetSceneCount()
         << ' ' << fName << " scene handlers extanct." << G4endl;
  return pScene;
}

G4VViewer* G4HepRepFile::CreateViewer(G4VSceneHandler& scene,
			       const G4String& name) {
  G4VViewer* pView =
    new G4HepRepFileViewer((G4HepRepFileSceneHandler&) scene, name);
  if (pView) {
    if (pView->GetViewId() < 0) {
      G4cout <<
	"G4HepRepFile::CreateViewer: ERROR flagged by negative"
        " view id in G4HepRepFileViewer creation."
        "\n Destroying view and returning null pointer."
	     << G4endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cout <<
      "G4HepRepFile::CreateViewer: ERROR: null pointer on new G4HepRepFileViewer."
	   << G4endl;
  }
  return pView;
}

HepRepXMLWriter* G4HepRepFile::GetHepRepXMLWriter () {
    return hepRepXMLWriter;
}
