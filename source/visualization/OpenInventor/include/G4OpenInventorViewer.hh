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

#ifndef G4OPENINVENTORVIEWER_HH
#define G4OPENINVENTORVIEWER_HH

#ifdef G4VIS_BUILD_OI_DRIVER

// Inheritance :
#include "G4VViewer.hh"

class SoSelection;
class Geant4_SoImageWriter;
class Geant4_SoGL2PSAction;

class G4OpenInventorSceneHandler;
class G4VInteractorManager;

//
// Base class for various OpenInventorView classes.
//
class G4OpenInventorViewer: public G4VViewer {
public: //G4VViewer
  virtual void DrawView();
  virtual void ShowView();
  virtual void ClearView();
  virtual void SetView();
  virtual void KernelVisitDecision();
public:
  G4OpenInventorViewer(G4OpenInventorSceneHandler& scene,
		       const G4String& name = "");
  virtual ~G4OpenInventorViewer();
protected:
  virtual void ViewerRender() = 0;
  void Escape();
  void WritePostScript(const G4String& file = "g4out.ps");
  void WritePixmapPostScript(const G4String& file = "g4out.ps");
  void WriteInventor(const G4String& file = "g4out.iv");
  void CountTriangles();
private:
  G4bool CompareForKernelVisit(G4ViewParameters&);
private:
  G4ViewParameters fLastVP;  // Memory for making kernel visit decisions.
protected:
  G4OpenInventorSceneHandler& fG4OpenInventorSceneHandler;
  G4VInteractorManager* fInteractorManager;
  SoSelection* fSoSelection;
  Geant4_SoImageWriter* fSoImageWriter;
  Geant4_SoGL2PSAction* fGL2PSAction;
};

#endif

#endif
