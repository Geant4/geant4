// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FukuiRendererView.hh,v 1.1 1999-01-07 16:14:37 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Satoshi TANAKA, Fri Jun 28 12:10:14 JST 1996
// FukuiRenderer view - opens window, hard copy, etc.

//=================//
#ifdef G4VIS_BUILD_DAWN_DRIVER
//=================//


#ifndef G4FUKUI_RENDERER_VIEW_HH
#define G4FUKUI_RENDERER_VIEW_HH

#include "G4VView.hh"
#include "globals.hh"

class G4FukuiRendererScene ;

class G4FukuiRendererView: public G4VView {

	enum FRDEV {FRDEV_PS=1, FRDEV_XWIN=2, FRDEV_PS2=3, FRDEV_XWIN2=4, FRDEV_OPEN_GL=5, FRDEV_DEVICE_END=6} ;

public:
	//----- constructor and destructor
  G4FukuiRendererView  (G4FukuiRendererScene& scene, const G4String& name);
  ~G4FukuiRendererView ();

	//----- overriding base class methods
  void SetView   (); // Do nothing. SendViewParameters will do its job. 
  void ClearView ();
  void DrawView  ();
  void ShowView  ();
  void FlushView (); // ShowView() without calling Wait()

	//---- methods inherent to this class
  void Wait() ;
  void SendViewParameters ()  ;

private:
  void  SendDevice( FRDEV dev );
  void  SendDrawingStyle() ; 
  G4FukuiRendererScene& fScene; // Reference to Graphics Scene for this view.
};

#endif
#endif //G4VIS_BUILD_DAWN_DRIVER
