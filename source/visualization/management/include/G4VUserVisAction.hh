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

// Class Description:
//
// G4VUserVisAction is brought into effect by the command
// /vis/scene/add/userAction. The vis manager instantiates a
// G4CallbackModel when the vis action is registered, and adds
// it to the scene. The pure virtual Draw() method
// is invoked whenever the viewer needs to "visit the kernel", e.g.,
// to remake its graphical database, if any, or simply to refresh the
// screen. It is not intended to be called directly. It is called
// via the operator() which is defined to satisfy the template
// G4CallbackModel<G4VUserVisAction>.
//
// A concrete Draw() method would normally make use of the Draw() methods
// of the vis manager, e.g:
//   G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
//   if (pVisManager) {
//     pVisManager->Draw(G4Box("box",2*cm,2*cm,2*cm),
//                       G4VisAttributes(G4Colour(1,1,0)));
//   ...
// but it can also use pointers fpSceneHandler and pMP that give it
// direct access to the current scene handler and modeling parameters
// for more advanced use.
//
// See the User Guide for Application Developers.

#ifndef G4VUSERVISACTION_HH
#define G4VUSERVISACTION_HH

class G4VGraphicsScene;
class G4ModelingParameters;

class G4VUserVisAction
{
public: // With description
  G4VUserVisAction(): fpSceneHandler(nullptr), fpMP(nullptr) {}
  virtual ~G4VUserVisAction() {}
  void operator()(G4VGraphicsScene& scene, const G4ModelingParameters* pMP) {
    fpSceneHandler = &scene;
    fpMP           = pMP;
    Draw();
  }
protected:
  virtual void Draw() = 0;
  G4VGraphicsScene* fpSceneHandler;
  const G4ModelingParameters* fpMP;
};

#endif
