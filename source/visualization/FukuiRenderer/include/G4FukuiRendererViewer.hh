// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FukuiRendererViewer.hh,v 1.5 1999-12-15 14:54:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Satoshi TANAKA, Fri Jun 28 12:10:14 JST 1996
// FukuiRenderer viewer - opens window, hard copy, etc.

//=================//
#ifdef G4VIS_BUILD_DAWN_DRIVER
//=================//


#ifndef G4FUKUI_RENDERER_VIEWER_HH
#define G4FUKUI_RENDERER_VIEWER_HH

#include "G4VViewer.hh"
#include "globals.hh"

class G4FukuiRendererSceneHandler ;


class G4FukuiRendererViewer: public G4VViewer {

public:

	enum FRDEV {FRDEV_PS=1, FRDEV_XWIN=2, FRDEV_PS2=3, FRDEV_XWIN2=4, FRDEV_OPEN_GL=5, FRDEV_DEVICE_END=6} ;

	//----- constructor and destructor
  G4FukuiRendererViewer  (G4FukuiRendererSceneHandler& scene, const G4String& name);
  virtual ~G4FukuiRendererViewer ();

	//----- overriding base class methods
  void SetView   (); // Do nothing. SendViewParameters will do its job. 
  void ClearView ();
  void DrawView  ();
  void ShowView  ();

	//---- methods inherent to this class
  void Wait() ;
  void SendViewParameters ()  ;

  void  SendDevice( FRDEV dev );
  void  SendDrawingStyle() ; 

private:

  G4FukuiRendererSceneHandler& fSceneHandler; // Reference to Graphics Scene for this view.
};

#endif
#endif //G4VIS_BUILD_DAWN_DRIVER
