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

/**
 * @author Mark Donszelmann
 * @version $Id: G4HepRep.hh,v 1.6 2002-11-13 19:04:31 duns Exp $
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
