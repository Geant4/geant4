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
// $Id: G4DAWNFILEViewer.hh,v 1.7 2001-07-11 10:08:44 gunter Exp $
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

	//---- methods inherent to this class
  void SendViewParameters ()  ;
  const char* GetG4PrimViewer          () { return fG4PrimViewer ;}
  const char* GetG4PrimViewerInvocation() { return fG4PrimViewerInvocation ;}
  const char* GetPSViewer              () { return fPSViewer ;}
  void  SendDrawingStyleToDAWNGUI( G4std::ostream& out ) ;

private:
  G4DAWNFILESceneHandler& fSceneHandler; // Reference to Graphics Scene for this view.

  char  fG4PrimViewer          [32] ;
  char  fG4PrimViewerInvocation[64] ;
  char  fPSViewer              [32] ;

};

#endif
#endif //G4VIS_BUILD_DAWNFILE_DRIVER
