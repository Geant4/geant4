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
// Laurent Garnier  27th October 2011

#include "G4OpenGLStoredQtSceneHandler.hh"

#include "G4PhysicalVolumeModel.hh"
#include "G4LogicalVolumeModel.hh"
#include "G4Text.hh"
#include "G4VPhysicalVolume.hh"
#include "G4OpenGLQtViewer.hh"
#include <typeinfo>
#include <sstream>

G4OpenGLStoredQtSceneHandler::G4OpenGLStoredQtSceneHandler
(G4VGraphicsSystem& system,
 const G4String& name):
G4OpenGLStoredSceneHandler (system, name)
{}

G4OpenGLStoredQtSceneHandler::~G4OpenGLStoredQtSceneHandler ()
{}

G4bool G4OpenGLStoredQtSceneHandler::ExtraPOProcessing
(const G4Visible& visible, size_t currentPOListIndex)
{
  G4bool usesGLCommands = true;

  try {
    const G4Text& g4Text = dynamic_cast<const G4Text&>(visible);
    G4TextPlus* pG4TextPlus = new G4TextPlus(g4Text);
    pG4TextPlus->fProcessing2D = fProcessing2D;
    fPOList[currentPOListIndex].fpG4TextPlus = pG4TextPlus;
    usesGLCommands = false;
  }
  catch (const std::bad_cast&) {}  // No special action if not text.  Just carry on.

  G4PhysicalVolumeModel* pPVModel =
    dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  G4LogicalVolumeModel* pLVModel =
    dynamic_cast<G4LogicalVolumeModel*>(pPVModel);
  if (pPVModel && !pLVModel) {

    // This call comes from a G4PhysicalVolumeModel.  drawnPVPath is
    // the path of the current drawn (non-culled) volume in terms of
    // drawn (non-culled) ancestors.  Each node is identified by a
    // PVNodeID object, which is a physical volume and copy number.  It
    // is a vector of PVNodeIDs corresponding to the geometry hierarchy
    // actually selected, i.e., not culled.
    //    typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
    //    typedef std::vector<PVNodeID> PVPath;

    // The simplest algorithm, used by the Open Inventor Driver
    // developers, is to rely on the fact the G4PhysicalVolumeModel
    // traverses the geometry hierarchy in an orderly manner.  The last
    // mother, if any, will be the node to which the volume should be
    // added.  So it is enough to keep a map of scene graph nodes keyed
    // on the volume path ID.  Actually, it is enough to use the logical
    // volume as the key.  (An alternative would be to keep the PVNodeID
    // in the tree and match the PVPath from the root down.)

    // BUT IN OPENGL, IF THERE ARE TRANSPARENT OBJECTS, VOLUMES DO NOT
    // ARRIVE IN THE ABOVE ORDER.  (TRANSPARENT OBJECTS ARE DRWAN
    // LAST.)  SO WE MUST BE MORE SOPHISTICATED IN CONSTRUCTING A
    // TREE.

    // build a path for tree viewer
    G4OpenGLQtViewer* pGLViewer = dynamic_cast<G4OpenGLQtViewer*>(fpViewer);
    if ( pGLViewer ) {
      pGLViewer->addPVSceneTreeElement(fpModel->GetCurrentDescription(),pPVModel,(G4int)currentPOListIndex);
    }

  } else {  // Not from a G4PhysicalVolumeModel.

    if (fpModel) {

      
      // build a path for tree viewer
      G4OpenGLQtViewer* pGLViewer = dynamic_cast<G4OpenGLQtViewer*>(fpViewer);
      if ( pGLViewer ) {
        pGLViewer->addNonPVSceneTreeElement(fpModel->GetType(),(G4int)currentPOListIndex,fpModel->GetCurrentDescription().data(),visible);
      }
    }
  }

  return usesGLCommands;
}

G4bool G4OpenGLStoredQtSceneHandler::ExtraTOProcessing
(const G4Visible& visible, size_t currentTOListIndex)
{

  G4bool usesGLCommands = true;

  try {
    const G4Text& g4Text = dynamic_cast<const G4Text&>(visible);
    G4TextPlus* pG4TextPlus = new G4TextPlus(g4Text);
    pG4TextPlus->fProcessing2D = fProcessing2D;
    fTOList[currentTOListIndex].fpG4TextPlus = pG4TextPlus;
    usesGLCommands = false;
  }
  catch (const std::bad_cast&) {}  // Do nothing if not text.

  return usesGLCommands;
}

void G4OpenGLStoredQtSceneHandler::ClearStore () {

  //G4cout << "G4OpenGLStoredQtSceneHandler::ClearStore" << G4endl;

  G4OpenGLStoredSceneHandler::ClearStore ();  // Sets need kernel visit, etc.
  // Should recreate the tree
  G4OpenGLQtViewer* pGLQtViewer = dynamic_cast<G4OpenGLQtViewer*>(fpViewer);
  if ( pGLQtViewer ) {
    pGLQtViewer->clearTreeWidget();
  }
}

void G4OpenGLStoredQtSceneHandler::ClearTransientStore () {

  //G4cout << "G4OpenGLStoredQtSceneHandler::ClearTransientStore" << G4endl;

  G4OpenGLStoredSceneHandler::ClearTransientStore ();

  // Should recreate the tree
  // Make sure screen corresponds to graphical database...
  // FIXME : L.Garnier April 2012 : Could cause a infinite loop ?
  if (fpViewer) {
    fpViewer -> SetView ();
    fpViewer -> ClearView ();
    fpViewer -> DrawView ();
  }
}

void G4OpenGLStoredQtSceneHandler::SetScene(G4Scene* pScene){

  if (pScene != fpScene) {
    G4OpenGLQtViewer* pGLQtViewer = dynamic_cast<G4OpenGLQtViewer*>(fpViewer);
    if ( pGLQtViewer ) {
      pGLQtViewer->clearTreeWidget();
    }
  }
  G4VSceneHandler::SetScene(pScene);
}
