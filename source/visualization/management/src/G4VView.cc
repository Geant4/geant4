// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VView.cc,v 1.2 1999-01-08 16:33:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  27th March 1996
// Abstract interface class for graphics views.

#include "G4VView.hh"

#include "G4ios.hh"
#ifdef WIN32
#include <strstrea.h>
#else
#include <strstream.h>
#endif

#include "G4VisManager.hh"
#include "G4VGraphicsSystem.hh"
#include "G4VScene.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4Event.hh"

G4VView::G4VView (G4VScene& scene, G4int id, const G4String& name):
fScene (scene),
fViewId (id),
fModified (true),
fNeedKernelVisit (true)
{
  G4VisManager* pVMan = G4VisManager::GetInstance ();
  fVP = pVMan -> GetCurrentViewParameters ();
  if (name == "") {
    char charname [50];
    ostrstream ost (charname, 50);
    ost << fScene.GetName () << '-' << fViewId << ends;
    fName = charname;
  }
  else {
    fName = name;
  }
}

G4VView::~G4VView () {}

const G4VisAttributes* G4VView::GetApplicableVisAttributes
(const G4VisAttributes* pVisAttribs) const {
  // If the pointer is null, pick up the default vis attributes from
  // the view parameters.
  if (!pVisAttribs)
    pVisAttribs = GetViewParameters ().GetDefaultVisAttributes ();
  return pVisAttribs;
}

void G4VView::NeedKernelVisit () {
  // Notify all views that a kernel visit is required.
  const G4VViewList& viewList = fScene.GetViewList ();
  for (int i = 0; i < viewList.entries (); i++) {
    viewList [i] -> SetNeedKernelVisit ();
  }
}

void G4VView::ShowView () {}

void G4VView::ProcessView () {

  // If view parameters have been modified, SetView () works out consequences. 
  if (fModified) {
    fModified = false;
    SetView ();
  }

  // If the scene data has changed (fNeedVisit is set to true in
  // G4VScene::SetSceneData), or if the concrete view has decided that it
  // necessary to visit the kernel (this should be Done in the concrete
  // object's DrawView ())...
  if (fNeedKernelVisit) {
    fNeedKernelVisit = false;
    fScene.ProcessScene (*this);
  }

  /*********************************************

Note that hits and digi drawing is totally the responsibility of the
"user" at present.  So the following code just represents an
alternative approach, namely, the user sets a flag, then the Vis
Manager draws them.  I will leave the code hear for now; also the
flags and methods in G4ViewParameters.

  // If event data is available...
  G4Event* pEvent = G4VisManager::GetInstance () -> GetEvent ();
  if (pEvent) {
    if (fVP.IsViewHits ()) {
      G4cout << "Drawing hits..." << endl;
    }
    if (fVP.IsViewDigis ()) {
      G4cout << "Drawing digis..." << endl;
    }
  }
  ***********************************************/
}

void G4VView::SetViewParameters (const G4ViewParameters& vp) {
  fVP = vp;
  fModified = true;
}

ostream& operator << (ostream& os, const G4VView& v) {
  os << "View " << v.fName << ":\n";
  os << v.fVP;
  return os;
}
