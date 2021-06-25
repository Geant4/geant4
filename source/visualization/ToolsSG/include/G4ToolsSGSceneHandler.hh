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
// John Allison  6th October 2020

#if defined (G4VIS_BUILD_TOOLSSG_DRIVER) || defined (G4VIS_USE_TOOLSSG)

#ifndef G4TOOLSSGSCENEHANDLER_HH
#define G4TOOLSSGSCENEHANDLER_HH

#include "G4VSceneHandler.hh"

#include <vector>

namespace tools {namespace sg {class separator;}}
class G4ToolsSGNode;

class G4ToolsSGSceneHandler: public G4VSceneHandler {
  typedef G4VSceneHandler parent;
public:  
  G4ToolsSGSceneHandler(G4VGraphicsSystem& system,const G4String& name);
  virtual ~G4ToolsSGSceneHandler();
protected:  
  //G4ToolsSGSceneHandler(const G4ToolsSGSceneHandler& a_from):parent(a_from){}
  G4ToolsSGSceneHandler& operator=(const G4ToolsSGSceneHandler&){return *this;}
public:
  
  void PreAddSolid(const G4Transform3D& objectTransformation,
                   const G4VisAttributes&);
  void PostAddSolid();

  void BeginPrimitives2D(const G4Transform3D& objectTransformation);
  void EndPrimitives2D();

  void BeginPrimitives(const G4Transform3D& objectTransformation);
  void EndPrimitives();

  virtual void AddPrimitive(const G4Polyline&);
  virtual void AddPrimitive(const G4Scale&);
  virtual void AddPrimitive(const G4Text&);
  virtual void AddPrimitive(const G4Circle&);
  virtual void AddPrimitive(const G4Square&);
  virtual void AddPrimitive(const G4Polymarker&);
  virtual void AddPrimitive(const G4Polyhedron&);

  void ClearStore ();
  void ClearTransientStore ();
  
public:
  tools::sg::separator* ToolsSGScene() {return fpToolsSGScene;}
protected:

  void EstablishBaseNodes();

  tools::sg::separator* GetOrCreateNode();  // For next solid or primitive

  static G4int fSceneIdCount;  // Counter for Qt3D scene handlers.

  tools::sg::separator* fpToolsSGScene;
  tools::sg::separator* fpTransientObjects;
  tools::sg::separator* fpPersistentObjects;
  std::vector<G4ToolsSGNode*> fpPhysicalVolumeObjects;  // Multiple worlds

};

#endif

#endif
