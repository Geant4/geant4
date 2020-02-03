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
/// \file visualization/perspective/include/PerspectiveVisAction.hh
/// \brief Definition of the PerspectiveVisAction class
//
//

#ifndef PERSPECTIVEVISACTION_HH
#define PERSPECTIVEVISACTION_HH

#include "G4VUserVisAction.hh"

#include "G4String.hh"
#include <map>
#include <vector>

/*
class G4AttDef;
class G4AttValue;
*/

class G4VVisManager;
class G4VSolid;
class G4VisAttributes;

class PerspectiveVisAction: public G4VUserVisAction {
public:
  PerspectiveVisAction();
  void SetOptionString(const G4String& optionString)
  {fOptionString = optionString;}
  void SetScene(const G4String& scene)
  {fScene = scene;}
  virtual void Draw();
private:
  void ExtendedDraw (const G4VSolid&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D());
  void RoomAndChair();
  void Chair(const G4VisAttributes&, const G4Transform3D&);
  G4VVisManager* fpVisManager;
  G4String fOptionString;
  G4String fScene;
  G4double fRoomX, fRoomY, fRoomZ,  // Half lengths.
    fWindowX, fWindowY, fWindowZ, fWindowSillHeight, fWindowOffset,
    fDoorFrameX, fDoorFrameY, fDoorFrameZ, fDoorFrameOffset,
    fDoorX, fDoorY, fDoorZ,
    fChairX,          // Half width.
    fChairY,          // Half depth.
    fChairZ,          // Half height.
    fChairSeat,       // Half height of top of seat.
    fChairThickness;  // Half thicknes of back, seat, legs.
};

#endif
