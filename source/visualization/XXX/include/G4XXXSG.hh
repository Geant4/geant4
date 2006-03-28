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
// $Id: G4XXXSG.hh,v 1.1 2006-03-28 17:16:41 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  10th March 2005
// A template for a sophisticated graphics driver with a scene graph.
//?? Lines or sections marked like this require specialisation for your driver.

#if defined (G4VIS_BUILD_XXXSG_DRIVER) || (G4VIS_USE_XXXSG)

#ifndef G4XXXSG_HH
#define G4XXXSG_HH

#include "G4VGraphicsSystem.hh"

class G4XXXSG: public G4VGraphicsSystem {
public:
  G4XXXSG();
  virtual ~G4XXXSG();
  G4VSceneHandler* CreateSceneHandler(const G4String& name = "");
  G4VViewer* CreateViewer (G4VSceneHandler&, const G4String& name = "");
};

#endif

#endif
