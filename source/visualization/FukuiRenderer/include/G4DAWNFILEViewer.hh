// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DAWNFILEViewer.hh,v 1.3 1999-05-10 15:38:28 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Satoshi TANAKA
// DAWNFILE viewer - opens window, hard copy, etc.

//=================//
#ifdef G4VIS_BUILD_DAWNFILE_DRIVER
//=================//


#ifndef G4DAWNFILE_VIEWER_HH
#define G4DAWNFILE_VIEWER_HH

#include "G4VViewer.hh"
#include "globals.hh"

class G4DAWNFILESceneHandler ;

class G4DAWNFILEViewer: public G4VViewer {

	enum FRDEV {FRDEV_PS=1, FRDEV_XWIN=2, FRDEV_PS2=3, FRDEV_XWIN2=4, FRDEV_OPEN_GL=5, FRDEV_DEVICE_END=6} ;

public:
	//----- constructor and destructor
  G4DAWNFILEViewer  (G4DAWNFILESceneHandler& scene, const G4String& name = "");
  virtual ~G4DAWNFILEViewer ();

	//----- overriding base class methods
  void SetView   (); // Do nothing. SendViewParameters will do its job. 
  void ClearView ();
  void DrawView  ();
  void ShowView  ();
  void FlushView (); // ShowView() without calling Wait()

	//---- methods inherent to this class
  void SendViewParameters ()  ;
  const char* GetG4PrimViewer          () { return fG4PrimViewer ;}
  const char* GetG4PrimViewerInvocation() { return fG4PrimViewerInvocation ;}
  const char* GetPSViewer              () { return fPSViewer ;}
  void  SendDrawingStyleToDAWNGUI( ostream& out ) ;

private:
  G4DAWNFILESceneHandler& fSceneHandler; // Reference to Graphics Scene for this view.

  char  fG4PrimViewer          [32] ;
  char  fG4PrimViewerInvocation[64] ;
  char  fPSViewer              [32] ;

};

#endif
#endif //G4VIS_BUILD_DAWNFILE_DRIVER
