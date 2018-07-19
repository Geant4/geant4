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
// 
// Jeff Kallenbach 01 Aug 1996
// OpenInventor viewer - opens window, hard copy, etc.

#ifndef G4OPENINVENTORXTEXTENDEDVIEWER_HH
#define G4OPENINVENTORXTEXTENDEDVIEWER_HH

#ifdef G4VIS_BUILD_OI_DRIVER

// Inheritance :
#include "G4OpenInventorViewer.hh"

#include <X11/Intrinsic.h>
#include <Inventor/nodes/SoEventCallback.h>

class G4OpenInventorXtExaminerViewer;

class G4OpenInventorXtExtendedViewer: public G4OpenInventorViewer {
public: //G4VViewer
  virtual void FinishView();
  virtual void SetView();
protected:
  virtual void ViewerRender();
  virtual SoCamera* GetCamera();
public:
  G4OpenInventorXtExtendedViewer(G4OpenInventorSceneHandler& scene,
		         const G4String& name = "");
  void Initialise();

  virtual ~G4OpenInventorXtExtendedViewer();
private:
  Widget AddMenu(Widget,const G4String&,const G4String&);
  void AddButton(Widget,const G4String&,XtCallbackProc);
private:
  static void PostScriptCbk(Widget,XtPointer,XtPointer);
  static void PDFCbk(Widget,XtPointer,XtPointer);
  static void PixmapPostScriptCbk(Widget,XtPointer,XtPointer);
  static void WriteInventorCbk(Widget,XtPointer,XtPointer);
  static void LoadInventorCbk(Widget,XtPointer,XtPointer);
  static void EscapeCbk(Widget,XtPointer,XtPointer);
  static void EscapeFromKeyboardCbk(void * o);
  static void SceneGraphStatisticsCbk(Widget,XtPointer,XtPointer);
  static void EraseDetectorCbk(Widget,XtPointer,XtPointer);
  static void EraseEventCbk(Widget,XtPointer,XtPointer);
  static void SetSolidCbk(Widget,XtPointer,XtPointer);
  static void SetWireFrameCbk(Widget,XtPointer,XtPointer);
  static void SetReducedWireFrameCbk(Widget,XtPointer,XtPointer);
  static void SetFullWireFrameCbk(Widget,XtPointer,XtPointer);
  static void UpdateSceneCbk(Widget,XtPointer,XtPointer);
  static void HelpCbk(Widget,XtPointer,XtPointer);
  static void HelpCancelCbk(Widget,XtPointer,XtPointer);
  static void SetPreviewCbk(Widget,XtPointer,XtPointer);
  static void SetPreviewAndFullCbk(Widget,XtPointer,XtPointer);
private:
  Widget fShell;
  G4OpenInventorXtExaminerViewer* fViewer;
  Widget fHelpForm;
  Widget fHelpText;
};

#endif

#endif
