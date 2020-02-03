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

// Frederick Jones TRIUMF 07 November 2017
// Open Inventor viewer using SoQt.

#ifndef G4OPENINVENTORQTVIEWER_HH
#define G4OPENINVENTORQTVIEWER_HH

#ifdef G4VIS_BUILD_OIQT_DRIVER

// Inheritance :
#include "G4OpenInventorViewer.hh"

// This causes clash with enum CursorShape in qnamespace.h !
//#include <X11/Intrinsic.h>

#include <Inventor/nodes/SoEventCallback.h>

//class SoQtExaminerViewer;
class G4OpenInventorQtExaminerViewer;

class G4OpenInventorQtViewer: public G4OpenInventorViewer {
public:
  G4OpenInventorQtViewer(G4OpenInventorSceneHandler& scene,
		         const G4String& name = "");
  virtual ~G4OpenInventorQtViewer();
  void Initialise();
public: //G4VViewer
  virtual void FinishView();
  virtual void SetView();
protected:
  virtual void ViewerRender();
  virtual SoCamera* GetCamera();

protected:

   //SoQtExaminerViewer* fViewer;
  G4OpenInventorQtExaminerViewer* fViewer;

  // FWJ from OpenInventorXtViewer.hh
  // Widget AddMenu(Widget,const G4String&,const G4String&);
  // void AddButton(Widget,const G4String&,XtCallbackProc);
  // static void PostScriptCbk(Widget,XtPointer,XtPointer);
  // static void PixmapPostScriptCbk(Widget,XtPointer,XtPointer);
  // static void WriteInventorCbk(Widget,XtPointer,XtPointer);
  // static void EscapeCbk(Widget,XtPointer,XtPointer);
  // static void SceneGraphStatisticsCbk(Widget,XtPointer,XtPointer);
  // static void EraseDetectorCbk(Widget,XtPointer,XtPointer);
  // static void EraseEventCbk(Widget,XtPointer,XtPointer);
  // static void SetSolidCbk(Widget,XtPointer,XtPointer);
  // static void SetWireFrameCbk(Widget,XtPointer,XtPointer);
  // static void SetReducedWireFrameCbk(Widget,XtPointer,XtPointer);
  // static void SetFullWireFrameCbk(Widget,XtPointer,XtPointer);
  // static void UpdateSceneCbk(Widget,XtPointer,XtPointer);
  // static void HelpCbk(Widget,XtPointer,XtPointer);
  // static void HelpCancelCbk(Widget,XtPointer,XtPointer);
  // static void SetPreviewCbk(Widget,XtPointer,XtPointer);
  // static void SetPreviewAndFullCbk(Widget,XtPointer,XtPointer);
  // Widget fShell;
  //  Widget fHelpForm;
  //  Widget fHelpForm;
  //  Widget fHelpText;
};

#endif

#endif
