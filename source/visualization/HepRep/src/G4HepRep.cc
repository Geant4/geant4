/**
 * @author Mark Donszelmann
 * @version $Id: G4HepRep.cc,v 1.4 2002-11-13 18:39:36 duns Exp $
 */

//G4
#include "G4Types.hh"
#include "G4HepRepViewer.hh"
#include "G4HepRepSceneHandler.hh"

//This
#include "G4HepRep.hh"

using namespace HEPREP;
using namespace std;

G4HepRep::G4HepRep ()
        : G4VGraphicsSystem ("HepRep",
        "HepRep Generic Driver for XML, RMI and CORBA",
        G4VGraphicsSystem::threeD) {
}

G4HepRep::~G4HepRep () {
}

G4VSceneHandler* G4HepRep::CreateSceneHandler (const G4String& name) {
    G4HepRepSceneHandler* scene = new G4HepRepSceneHandler (*this, name);
    return scene;
}

G4VViewer* G4HepRep::CreateViewer (G4VSceneHandler& scene, const G4String& name) {
    G4VViewer* view  = new G4HepRepViewer ((G4HepRepSceneHandler&)scene, name);
    if (view) {
        if (view->GetViewId() < 0) {
            G4cout << "G4HepRep::CreateViewer: ERROR flagged by negative view ID." << G4endl;
            delete view;
            view = NULL;
        }
    } else {
        G4cout << "G4HepRep::CreateViewer: null pointer on new G4HepRepViewer." << G4endl;
    }
    return view;
}
