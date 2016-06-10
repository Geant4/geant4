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
// $Id: G4HepRepFile.cc 78838 2014-01-28 08:46:17Z gcosmo $
//
// 
// Joseph Perl  1 October 2001
// A graphics system to dump geometry hierarchy to the
// HepRep graphics format (HepRep version 1).

//HepRep
#include "G4HepRepFileXMLWriter.hh"

#include "G4HepRepFile.hh"
#include "G4HepRepFileSceneHandler.hh"
#include "G4HepRepFileViewer.hh"
#include "G4HepRepMessenger.hh"

static G4HepRepFileXMLWriter* hepRepXMLWriter;

G4HepRepFile::G4HepRepFile():
  G4VGraphicsSystem("G4HepRepFile",
		    "HepRepFile",
		    "A HepRep (format 1) ascii file driver",
		    G4VGraphicsSystem::fileWriter) {
		G4HepRepMessenger::GetInstance();
        hepRepXMLWriter = new G4HepRepFileXMLWriter();
}

G4HepRepFile::~G4HepRepFile()
{
        delete hepRepXMLWriter;
}

G4VSceneHandler* G4HepRepFile::CreateSceneHandler(const G4String& name) {
  G4VSceneHandler* pScene = new G4HepRepFileSceneHandler(*this, name);
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

G4HepRepFileXMLWriter* G4HepRepFile::GetHepRepXMLWriter () {
    return hepRepXMLWriter;
}
