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
// $Id: G4XXX.hh,v 1.1 2001-08-17 23:04:25 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A graphics system to dump geometry hierarchy.

#ifndef G4XXX_HH
#define G4XXX_HH

#include "G4VGraphicsSystem.hh"

class G4XXX: public G4VGraphicsSystem {
public:
  G4XXX();
  virtual ~G4XXX();
  G4VSceneHandler* CreateSceneHandler(const G4String& name = "");
  G4VViewer* CreateViewer (G4VSceneHandler&, const G4String& name = "");
};

#endif
