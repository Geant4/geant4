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
// $Id: G4HepRep.cc,v 1.1 2001-08-24 23:06:37 perl Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A graphics system to dump geometry hierarchy.

//HepRep
#include "JHepRep.hh"
#include "JHepRepFactory.hh"

#include "G4HepRep.hh"
#include "G4HepRepSceneHandler.hh"
#include "G4HepRepViewer.hh"

static JHepRepFactory* factory;
static JHepRep* heprep = NULL;

G4HepRep::G4HepRep():
  G4VGraphicsSystem("G4HepRep",
		    "HepRep",
		    "A template graphics driver",
		    G4VGraphicsSystem::noFunctionality) {
        factory = new JHepRepFactory();
	heprep = factory->CreateHepRep();
	//factory->RunWired(heprep);
}

G4HepRep::~G4HepRep() {}

G4VSceneHandler* G4HepRep::CreateSceneHandler(const G4String& name) {
  G4VSceneHandler* pScene = new G4HepRepSceneHandler(*this, name);
  G4cout << G4HepRepSceneHandler::GetSceneCount()
         << ' ' << fName << " scene handlers extanct." << G4endl;
  return pScene;
}

G4VViewer* G4HepRep::CreateViewer(G4VSceneHandler& scene,
			       const G4String& name) {
  G4VViewer* pView =
    new G4HepRepViewer((G4HepRepSceneHandler&) scene, name);
  if (pView) {
    if (pView->GetViewId() < 0) {
      G4cout <<
	"G4HepRep::CreateViewer: ERROR flagged by negative"
        " view id in G4HepRepViewer creation."
        "\n Destroying view and returning null pointer."
	     << G4endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cout <<
      "G4HepRep::CreateViewer: ERROR: null pointer on new G4HepRepViewer."
	   << G4endl;
  }
  return pView;
}

JHepRepFactory* G4HepRep::GetHepRepFactory () {
    return factory;
}

JHepRep* G4HepRep::GetHepRep () {
    return heprep;
}
