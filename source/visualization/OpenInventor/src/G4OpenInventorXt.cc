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
// $Id: G4OpenInventorXt.cc,v 1.2 2004/11/22 15:00:48 gbarrand Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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

G4OpenInventorXt::G4OpenInventorXt ()
:G4OpenInventor("OpenInventorXt","OIX",G4VGraphicsSystem::threeD)
{
  SetInteractorManager (G4Xt::getInstance ());
  GetInteractorManager () -> 
    RemoveDispatcher((G4DispatchFunction)XtDispatchEvent);  
  GetInteractorManager () -> 
    AddDispatcher   ((G4DispatchFunction)SoXt::dispatchEvent);

  Widget top = (Widget)GetInteractorManager()->GetMainInteractor();


  if(getenv("XENVIRONMENT")==NULL) {
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

  SoXt::init(top);

  InitNodes();
}

G4OpenInventorXt::~G4OpenInventorXt () {}
G4VViewer* G4OpenInventorXt::CreateViewer (G4VSceneHandler& scene, const G4String& name) 
{
  G4OpenInventorSceneHandler* pScene = (G4OpenInventorSceneHandler*)&scene;
  return new G4OpenInventorXtViewer (*pScene, name);
}


#endif
