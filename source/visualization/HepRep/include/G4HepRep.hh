/**
 * @author Mark Donszelmann
 * @version $Id: G4HepRep.hh,v 1.4 2002-11-13 18:38:40 duns Exp $
 */
#ifndef G4HEPREP_HH
#define G4HEPREP_HH 1

//G4
#include "G4VGraphicsSystem.hh"

class G4HepRep: public G4VGraphicsSystem {
    public:
        G4HepRep ();
        virtual ~G4HepRep ();
        G4VSceneHandler* CreateSceneHandler (const G4String& name = "");
        G4VViewer* CreateViewer  (G4VSceneHandler&, const G4String& name = "");

};

#endif
