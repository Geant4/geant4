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

#ifdef G4VIS_BUILD_OIQT_DRIVER

#include "G4OpenInventorQtExaminerViewer.hh"

#include "G4ios.hh"

#include <Inventor/Qt/SoQt.h>
#include <Inventor/events/SoKeyboardEvent.h>


G4OpenInventorQtExaminerViewer* G4OpenInventorQtExaminerViewer::viewer = 0;


// Constructor
G4OpenInventorQtExaminerViewer::
G4OpenInventorQtExaminerViewer(QWidget* parent, const char* name, SbBool embed,
                               SoQtFullViewer::BuildFlag flag, 
                               SoQtViewer::Type type)
   : SoQtExaminerViewer(parent, name, embed, flag, type)
{
   G4cout << "G4OpenInventorQtExaminerViewer CONSTRUCTOR CALLED" << G4endl;
   viewer = this;
}

// Destructor
G4OpenInventorQtExaminerViewer::~G4OpenInventorQtExaminerViewer()
{
   //   if (superimposition != NULL) {
   //      removeSuperimposition(superimposition);
   //      superimposition->unref();
   //      superimposition = NULL;
   //   }
   //   if (animateSensor->isScheduled())
   //      animateSensor->unschedule();
   //   delete animateSensor;
   //   delete sceneChangeSensor;
   //
   //   delete[] curViewPtName;
   //   delete searcher;

   viewer = 0;
}


// FWJ try to handle events as in XtExaminerViewer
SbBool 
G4OpenInventorQtExaminerViewer::processSoEvent(const SoEvent* const ev)
{

   //   SoCamera *cam = getCamera();
   const SoType type(ev->getTypeId());

   if (type.isDerivedFrom(SoKeyboardEvent::getClassTypeId())) {
      SoKeyboardEvent* ke = (SoKeyboardEvent*)ev;

      if (SoKeyboardEvent::isKeyPressEvent(ev, ke->getKey())) {
         switch (ke->getKey()) {
         case SoKeyboardEvent::E:
            G4cout << "E KEY PRESSED, EXITING Qt MAIN LOOP" << G4endl;
            SoQt::exitMainLoop();
            break;
         default:
            break;
         }
      }
   }
   // Pass the event on to the viewer
   // Need some checks here as in Xt viewer?
   return SoQtExaminerViewer::processSoEvent(ev);
}

#endif
