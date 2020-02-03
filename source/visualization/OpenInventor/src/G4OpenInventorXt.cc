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
// OpenInventor graphics system factory.

#ifdef G4VIS_BUILD_OIX_DRIVER

// this :
#include "G4OpenInventorXt.hh"

#include <Inventor/Xt/SoXt.h>

#include "G4Xt.hh"
#include "G4OpenInventorSceneHandler.hh"
#include "G4OpenInventorXtViewer.hh"

// Work around for gcc8 Coverity cast warning
inline bool soxt_dispatch_event(void* a_event) {
   return SoXt::dispatchEvent((XEvent*)a_event);
}

G4OpenInventorXt::G4OpenInventorXt ()
:G4OpenInventor("OpenInventorXt","OIX",G4VGraphicsSystem::threeD)
,fInited(false)
{
}

void G4OpenInventorXt::Initialize()
{
  if(fInited) return; //Done

  SetInteractorManager (G4Xt::getInstance ());
  GetInteractorManager () -> 
     RemoveDispatcher(G4Xt::xt_dispatch_event);
  // Coverity gcc8 cast warning
  //    RemoveDispatcher((G4DispatchFunction)XtDispatchEvent);  
  GetInteractorManager () -> 
     AddDispatcher(soxt_dispatch_event);
  // Coverity gcc8 cast warning
  //    AddDispatcher   ((G4DispatchFunction)SoXt::dispatchEvent);

  Widget top = (Widget)GetInteractorManager()->GetMainInteractor();

  if(std::getenv("XENVIRONMENT")==NULL) {
    XrmDatabase database = XrmGetDatabase(XtDisplay(top));
    if(database!=NULL) {
      XrmPutLineResource(&database,"*topShadowColor:white");
      XrmPutLineResource(&database,"*bottomShadowColor:black");
      XrmPutLineResource(&database,"*foreground:black");
      XrmPutLineResource(&database,"*background:lightgrey");
      XrmPutLineResource(&database,"*borderColor:lightgrey");
      XrmPutLineResource(&database,"*fontList:-*-helvetica-bold-r-*-*-*-120-*-*-*-*-iso8859-1");
      XrmPutLineResource(&database,"*help_popup.title:Help");
      XrmPutLineResource(&database,"*helpCancel.labelString:Cancel");
      XrmPutLineResource(&database,"*helpText.editMode:multi_line_edit");
      XrmPutLineResource(&database,"*helpText.columns:60");
      XrmPutLineResource(&database,"*helpText.rows:20");
      XrmPutLineResource(&database,"*helpText.background:white");
      XrmPutLineResource(&database,"*helpText.fontList:*courier*-r-*--14-*");
      XrmPutLineResource(&database,"*helpText.maxLength:8000");
    }
  }

  if(!SoXt::getTopLevelWidget()) SoXt::init(top);

  InitNodes();

  fInited = true;
}

G4OpenInventorXt::~G4OpenInventorXt () {}
G4VViewer* G4OpenInventorXt::CreateViewer (G4VSceneHandler& scene, const G4String& name) 
{
  Initialize();
  G4OpenInventorSceneHandler* pScene = (G4OpenInventorSceneHandler*)&scene;
  return new G4OpenInventorXtViewer (*pScene, name);
}


#endif
