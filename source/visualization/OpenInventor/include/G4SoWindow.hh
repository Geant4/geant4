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
