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
// John Allison  17th June 2019

#if defined (G4VIS_BUILD_QT3D_DRIVER) || defined (G4VIS_USE_QT3D)

#ifndef G4QT3DVIEWER_HH
#define G4QT3DVIEWER_HH

#include "G4Qt3DSceneHandler.hh"

#include "G4VViewer.hh"

#include <Qt3DExtras>

class G4Qt3DViewer: public G4VViewer, public Qt3DExtras::Qt3DWindow
{
public:

  G4Qt3DViewer(G4Qt3DSceneHandler&,const G4String& name);
  virtual ~G4Qt3DViewer();
  void Initialise();
  void SetView();
  void ClearView();
  void DrawView();
  void ShowView();
  void FinishView();

  void SwitchToVisSubThread();
  void SwitchToMasterThread();

protected:

  void KernelVisitDecision ();
  G4bool CompareForKernelVisit(G4ViewParameters&);

  void keyPressEvent        (QKeyEvent*);
  void keyReleaseEvent      (QKeyEvent*);
  void mouseDoubleClickEvent(QMouseEvent*);
  void mouseMoveEvent       (QMouseEvent*);
  void mousePressEvent      (QMouseEvent*);
  void mouseReleaseEvent    (QMouseEvent*);
  void wheelEvent           (QWheelEvent*);

  G4ViewParameters fLastVP;  // Memory for making kernel visit decisions.
  G4Qt3DSceneHandler& fQt3DSceneHandler;

  QWidget* fUIWidget;

  G4bool fKeyPressed;
  int fKey;
  G4bool fMousePressed;
  G4double fMousePressedX, fMousePressedY;
};

#endif

#endif  // #if defined (G4VIS_BUILD_QT3D_DRIVER) || defined (G4VIS_USE_QT3D)
