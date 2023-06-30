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

#ifndef G4QT3DSCENEHANDLER_HH
#define G4QT3DSCENEHANDLER_HH

#include "G4VSceneHandler.hh"

#include <vector>

namespace Qt3DCore {
class QEntity;
}
class G4Qt3DQEntity;

class G4Qt3DSceneHandler: public G4VSceneHandler {
  
  friend class G4Qt3DViewer;

public:
  
  G4Qt3DSceneHandler(G4VGraphicsSystem& system,
                    const G4String& name);
  virtual ~G4Qt3DSceneHandler();
  
  void PreAddSolid(const G4Transform3D& objectTransformation,
                   const G4VisAttributes&);
  void PostAddSolid();

  void BeginPrimitives2D(const G4Transform3D& objectTransformation);
  void EndPrimitives2D();

  void BeginPrimitives(const G4Transform3D& objectTransformation);
  void EndPrimitives();

  using G4VSceneHandler::AddPrimitive;
  void AddPrimitive(const G4Polyline&);
  void AddPrimitive(const G4Polymarker&);
  void AddPrimitive(const G4Text&);
  void AddPrimitive(const G4Circle&);
  void AddPrimitive(const G4Square&);
  void AddPrimitive(const G4Polyhedron&);

  using G4VSceneHandler::AddCompound;
  void AddCompound(const G4Mesh&);

  void ClearStore ();
  void ClearTransientStore ();
  
protected:

  void EstablishG4Qt3DQEntities();
  G4Qt3DQEntity* CreateNewNode();  // For next solid or primitive

  static G4int fSceneIdCount;  // Counter for Qt3D scene handlers.

  Qt3DCore::QEntity* fpQt3DScene;
  Qt3DCore::QEntity* fpTransientObjects;
  Qt3DCore::QEntity* fpPersistentObjects;
  std::vector<G4Qt3DQEntity*> fpPhysicalVolumeObjects;  // Multiple worlds

};

#endif
