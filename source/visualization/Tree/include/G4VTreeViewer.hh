// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VTreeViewer.hh,v 1.3 2001-06-15 07:22:55 stanaka Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A base class for a dummy viewer to dump geometry hierarchy.

#ifndef G4VTREEVIEWER_HH
#define G4VTREEVIEWER_HH

#include "G4VViewer.hh"

class G4VTreeViewer: public G4VViewer {
public:
  G4VTreeViewer(G4VSceneHandler&,const G4String& name);
  virtual ~G4VTreeViewer();
  void SetView();
  void ClearView();
  void DrawView();
};

#endif
