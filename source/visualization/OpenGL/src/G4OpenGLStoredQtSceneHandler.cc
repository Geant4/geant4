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
// $Id: G4OpenGLStoredSceneHandler.cc,v 1.46 2010-11-10 17:11:20 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Laurent Garnier  27th October 2011

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#include "G4OpenGLStoredQtSceneHandler.hh"

#include "G4PhysicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4OpenGLQtViewer.hh"

G4OpenGLStoredQtSceneHandler::G4OpenGLStoredQtSceneHandler
(G4VGraphicsSystem& system,
 const G4String& name):
G4OpenGLStoredSceneHandler (system, name)
{}

G4OpenGLStoredQtSceneHandler::~G4OpenGLStoredQtSceneHandler ()
{}

void G4OpenGLStoredQtSceneHandler::ExtraPOProcessing
(size_t currentPOListIndex)
{

  G4PhysicalVolumeModel* pPVModel =
    dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  
  if (pPVModel) {

    // This call comes from a G4PhysicalVolumeModel.  drawnPVPath is
    // the path of the current drawn (non-culled) volume in terms of
    // drawn (non-culled) ancesters.  Each node is identified by a
    // PVNodeID object, which is a physical volume and copy number.  It
    // is a vector of PVNodeIDs corresponding to the geometry hierarchy
    // actually selected, i.e., not culled.
    typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
    typedef std::vector<PVNodeID> PVPath;
    const PVPath& drawnPVPath = pPVModel->GetDrawnPVPath();

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

    
    // Name, currentPOindex, copyNo
    std::vector < std::pair<std::string,std::pair <unsigned int, unsigned int> > > treeVect;
    std::pair<std::string,std::pair <unsigned int, unsigned int> > treeInfos;

    for (PVPath::const_iterator i = drawnPVPath.begin();
	 i != drawnPVPath.end(); ++i) {
      treeInfos.second.first = currentPOListIndex;
      treeInfos.second.second = i->GetCopyNo();
      treeInfos.first = i->GetPhysicalVolume()->GetName().data();
      treeVect.push_back(treeInfos);
    }

    // build a path for tree viewer
    G4OpenGLQtViewer* pGLViewer = dynamic_cast<G4OpenGLQtViewer*>(fpViewer);
    if ( pGLViewer ) {
      pGLViewer->addTreeElement(fpModel->GetCurrentDescription(),treeVect);
    }

  } else {  // Not from a G4PhysicalVolumeModel.

    if (fpModel) {
      // Create a place for current solid in root of scene graph tree...
      // Name, currentPOindex, copyNo
      std::vector < std::pair<std::string,std::pair <unsigned int, unsigned int> > > treeVect;
      std::pair<std::string,std::pair <unsigned int, unsigned int> > treeInfos;
      
      treeInfos.second.first = currentPOListIndex;
      treeInfos.second.second = 0;
      treeInfos.first = fpModel->GetCurrentTag().data();
      treeVect.push_back(treeInfos);
      
      // build a path for tree viewer
      G4OpenGLQtViewer* pGLViewer = dynamic_cast<G4OpenGLQtViewer*>(fpViewer);
      if ( pGLViewer ) {
        pGLViewer->addTreeElement(fpModel->GetCurrentDescription(),treeVect);
      }
    }
  }
}

void G4OpenGLStoredQtSceneHandler::ExtraTOProcessing
(size_t)
{
  //G4cout << "G4OpenGLStoredQtSceneHandler::ExtraTOProcessing: index: "
  //	 << currentTOListIndex << G4endl;
}

void G4OpenGLStoredQtSceneHandler::ClearStore () {

  //G4cout << "G4OpenGLStoredQtSceneHandler::ClearStore" << G4endl;

  G4OpenGLStoredSceneHandler::ClearStore ();  // Sets need kernel visit, etc.

  // Delete Qt Tree.
}

void G4OpenGLStoredQtSceneHandler::ClearTransientStore () {

  //G4cout << "G4OpenGLStoredQtSceneHandler::ClearTransientStore" << G4endl;

  G4OpenGLStoredSceneHandler::ClearTransientStore ();

  // Make sure screen corresponds to graphical database...
  if (fpViewer) {
    fpViewer -> SetView ();
    fpViewer -> ClearView ();
    fpViewer -> DrawView ();
  }
}

#endif
