//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
// Satoshi TANAKA
// DAWNFILE viewer - opens window, hard copy, etc.

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
  void  SendDrawingStyleToDAWNGUI( std::ostream& out ) ;

private:
  G4DAWNFILESceneHandler& fSceneHandler; // Reference to Graphics Scene for this view.

  char  fG4PrimViewer          [32] ;
  char  fG4PrimViewerInvocation[64] ;
  char  fPSViewer              [32] ;

};

#endif
