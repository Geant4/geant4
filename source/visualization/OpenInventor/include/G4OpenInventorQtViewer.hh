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

#if defined (G4VIS_BUILD_OIQT_DRIVER) || defined (G4VIS_USE_OIQT)

// Inheritance :
#include "G4OpenInventorViewer.hh"

#include <Inventor/nodes/SoEventCallback.h>

//class SoQtExaminerViewer;
class G4OpenInventorQtExaminerViewer;

#include <qobject.h>

class QMenuBar;
class QFont;
class QAction;


//class G4OpenInventorQtViewer: public G4OpenInventorViewer {
class G4OpenInventorQtViewer: public QObject,
                              public G4OpenInventorViewer {

  Q_OBJECT 

private Q_SLOTS :

   // File menu
   void FileWritePSCB();
   void FileWritePDFCB();
   void FileWriteIVCB();
   void FileEscapeCB();

   // Etc menu
   void EtcEraseDetectorCB();
   void EtcEraseEventCB();
   void EtcSetSolidCB();
   void EtcSetReducedWireframeCB();
   void EtcSetFullWireframeCB();
   void EtcVisibMInvisibDCB();
   void EtcVisibMVisibDCB();
   void EtcUpdateSceneCB();
   void EtcSceneGraphStatsCB();

private:

   // File menu
   QAction* FileWritePS;
   QAction* FileWritePDF;
   QAction* FileWriteIV;
   QAction* FileEscape;

   // Etc menu
   QAction* EtcEraseDetector;
   QAction* EtcEraseEvent;
   QAction* EtcSetSolid;
   QAction* EtcSetReducedWireframe;
   QAction* EtcSetFullWireframe;
   QAction* EtcVisibMInvisibD;
   QAction* EtcVisibMVisibD;
   QAction* EtcUpdateScene;
   QAction* EtcSceneGraphStats;

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

};

#endif

#endif
