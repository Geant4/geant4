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
// $Id: G4SceneHandlerList.hh,v 1.6 2001-07-11 10:09:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  May 1996

#ifndef G4SCENEHANDLERLIST_HH
#define G4SCENEHANDLERLIST_HH

#include "g4std/vector"
#include "G4VSceneHandler.hh"

class G4SceneHandlerList: public G4std::vector<G4VSceneHandler*> {
public:
  void remove(G4VSceneHandler*);
};

typedef G4std::vector<G4VSceneHandler*>::iterator G4SceneHandlerListIterator;
typedef G4std::vector<G4VSceneHandler*>::const_iterator
        G4SceneHandlerListConstIterator;

#endif
