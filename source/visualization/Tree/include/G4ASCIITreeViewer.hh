// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ASCIITreeViewer.hh,v 1.3 2001-06-15 07:22:55 stanaka Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A dummy viewer for ASCIITreeSceneHandler.

#ifndef G4ASCIITREEVIEWER_HH
#define G4ASCIITREEVIEWER_HH

#include "G4VTreeViewer.hh"

class G4ASCIITreeViewer: public G4VTreeViewer {
public:
  G4ASCIITreeViewer(G4VSceneHandler&,const G4String& name);
  virtual ~G4ASCIITreeViewer();
};

#endif
