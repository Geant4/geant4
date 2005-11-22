//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: PerspectiveVisAction.hh,v 1.1 2005-11-22 15:51:22 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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
  void Draw();
private:
  void ExtendedDraw (const G4VSolid&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D());
  void RoomAndChair();
  void Chair(const G4VisAttributes&, const G4Transform3D&);
  G4VVisManager* fpVisManager;
  G4String fOptionString;
  G4String fScene;
  G4double fRoomX, fRoomY, fRoomZ,  // Half lengths.
  fChairX,          // Half width.
  fChairY,          // Half depth.
  fChairZ,          // Half height.
  fChairSeat,       // Half height of top of seat.
  fChairThickness;  // Half thicknes of back, seat, legs.
};

#endif
