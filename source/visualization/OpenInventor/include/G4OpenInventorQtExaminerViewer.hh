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

// Frederick Jones TRIUMF 07 January 2018


#ifndef G4OPENINVENTORQTEXAMINERVIEWER_HH
#define G4OPENINVENTORQTEXAMINERVIEWER_HH

#ifdef G4VIS_BUILD_OIQT_DRIVER

//#include <map>
//#include <vector>
//#include <fstream>
//#include <Inventor/SbLinear.h>
//#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/events/SoKeyboardEvent.h>

class SoCoordinate3;
class SoFont;
class SoText2;
class SoPointSet;

class G4OpenInventorQtExaminerViewer : public SoQtExaminerViewer {

  //  friend class G4OpenInventorQtExaminerViewerMessenger;
  // FWJ
  friend class G4OpenInventorQtViewer;

private:

  static G4OpenInventorQtExaminerViewer* viewer;
   //  void (*escapeCallback)(void *);
   //  void * examinerObject;
   //  SbBool lshiftdown, rshiftdown, lctrldown, rctrldown;
  
public:

  // Same constructor as the SoQtExaminerViewer
  G4OpenInventorQtExaminerViewer(QWidget* parent = NULL,
	     const char* name = NULL,
	     SbBool embed = TRUE, 
	     SoQtFullViewer::BuildFlag flag = BUILD_ALL,
	     SoQtViewer::Type type = BROWSER);

  ~G4OpenInventorQtExaminerViewer();

protected:
  // FWJ why this?
  // Same constructor as the SoQtExaminerViewer 
  //  G4OpenInventorQtExaminerViewer(QWidget parent,
  //	     const char *name,
  //	     SbBool embed,
  //	     SoQtFullViewer::BuildFlag flag,
  //	     SoQtViewer::Type type,
  //	     SbBool build);

  // FWJ try to handle events as in XtExaminerViewer
  SbBool processSoEvent(const SoEvent* const event);

};

#endif

#endif /* G4OPENINVENTORQTEXAMINERVIEWER_HH */
