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
// $Id: G4OpenInventorXtViewer.hh,v 1.8 2004-11-22 11:09:58 gbarrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Jeff Kallenbach 01 Aug 1996
// OpenInventor viewer - opens window, hard copy, etc.

#ifndef G4OPENINVENTORXTVIEWER_HH
#define G4OPENINVENTORXTVIEWER_HH

#ifdef G4VIS_BUILD_OI_DRIVER

// Inheritance :
#include "G4OpenInventorViewer.hh"

#include <X11/Intrinsic.h>
class SoXtExaminerViewer;

class G4OpenInventorXtViewer: public G4OpenInventorViewer {
public: //G4VViewer
  virtual void FinishView();
protected:
  virtual void ViewerRender();
public:
  G4OpenInventorXtViewer(G4OpenInventorSceneHandler& scene,
		         const G4String& name = "");
  virtual ~G4OpenInventorXtViewer();
private:
  Widget AddMenu(Widget,const G4String&,const G4String&);
  void AddButton(Widget,const G4String&,XtCallbackProc);
private:
  static void PostScriptCbk(Widget,XtPointer,XtPointer);
  static void PixmapPostScriptCbk(Widget,XtPointer,XtPointer);
  static void WriteInventorCbk(Widget,XtPointer,XtPointer);
  static void EscapeCbk(Widget,XtPointer,XtPointer);
  static void SceneGraphStatisticsCbk(Widget,XtPointer,XtPointer);
  static void EraseDetectorCbk(Widget,XtPointer,XtPointer);
  static void EraseEventCbk(Widget,XtPointer,XtPointer);
  static void SetSolidCbk(Widget,XtPointer,XtPointer);
  static void SetWireFrameCbk(Widget,XtPointer,XtPointer);
  static void UpdateSceneCbk(Widget,XtPointer,XtPointer);
  static void HelpCbk(Widget,XtPointer,XtPointer);
private:
  Widget fShell;
  SoXtExaminerViewer* fViewer;
};

#endif

#endif
