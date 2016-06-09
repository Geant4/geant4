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
// $Id: G4OpenInventorWin.hh,v 1.1 2004/04/08 09:37:30 gbarrand Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// OpenInventor graphics system factory.

#if defined (G4VIS_BUILD_OIWIN32_DRIVER) || defined (G4VIS_USE_OIWIN32)

#ifndef G4OPENINVENTORWIN_HH
#define G4OPENINVENTORWIN_HH

#include "G4OpenInventor.hh"

class G4OpenInventorWin: public G4OpenInventor {
public:
  G4OpenInventorWin ();
  virtual ~G4OpenInventorWin ();
  G4VViewer* CreateViewer(G4VSceneHandler&,const G4String& name = "");
};

#endif

#endif
