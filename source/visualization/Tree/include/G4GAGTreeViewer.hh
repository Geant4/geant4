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
// Satoshi Tanaka  31th May 2001
// A dummy viewer for GAGTreeSceneHandler.

#ifndef G4GAGTREEVIEWER_HH
#define G4GAGTREEVIEWER_HH

#include "G4VTreeViewer.hh"

class G4GAGTreeViewer: public G4VTreeViewer {
public:
  G4GAGTreeViewer(G4VSceneHandler&,const G4String& name);
  virtual ~G4GAGTreeViewer();
};

#endif
