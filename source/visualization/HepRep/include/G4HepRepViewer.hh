/**
 * @author Mark Donszelmann
 * @version $Id: G4HepRepViewer.hh,v 1.4 2002-11-13 18:38:46 duns Exp $
 */

#ifndef G4HEPREPVIEWER_HH
#define G4HEPREPVIEWER_HH 1

// HepRep
#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepWriter.h"

// Geant4
#include "G4VViewer.hh"

//class G4HepRepSceneHandler;

class G4HepRepViewer: public G4VViewer {
    public:
        G4HepRepViewer (G4VSceneHandler& scene, const G4String& name = "");
        virtual ~G4HepRepViewer ();
        void SetView ();
        void ClearView ();
        void DrawView ();
        void ShowView ();
        void FinishView ();

    private:
        bool drawn;
};

#endif
