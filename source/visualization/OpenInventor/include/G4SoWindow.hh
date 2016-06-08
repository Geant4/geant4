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
#ifndef G4SOWINDOW_h
#define G4SOWINDOW_h 

#include <Inventor/Xt/SoXt.h>
#include <HEPVis/SbBasic.h>

class SbDict;

class G4SoWindow {
public:
  G4SoWindow(const char*,SbBool=FALSE);
  virtual ~G4SoWindow();
  Widget getWidget() const;
  void show();
  void hide();
  void setTitle(const char*);
  void setSize(const SbVec2s&);
  SbVec2s getSize();
  void setWMClose(SbBool);
  SbBool getWMClose() const;
  static void clearClass();
protected:
  void registerWidget(Widget);
  void unregisterWidget(Widget);
  static G4SoWindow* getWindow(Widget);
private:
  Widget fWidget;
  SbBool fWMClose;
  static SbDict* fWidgetDictionary;
#ifdef HEPVisXt
  static void widgetDestroyedCB(Widget,XtPointer,XtPointer);
public:
  static void closeWindow(Widget,XEvent*,char**,Cardinal*);
#endif
};

#endif
