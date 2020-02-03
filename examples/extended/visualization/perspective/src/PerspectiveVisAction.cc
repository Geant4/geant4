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
/// \file visualization/perspective/src/PerspectiveVisAction.cc
/// \brief Implementation of the PerspectiveVisAction class
//
//

#include "PerspectiveVisAction.hh"

#include "PerspectiveVisActionMessenger.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4Polyline.hh"
#include "G4Polyhedron.hh"
#include "G4Vector3D.hh"
#include "G4Point3D.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PerspectiveVisAction::PerspectiveVisAction():
  G4VUserVisAction(),
  fpVisManager(0),
  fOptionString("none"),
  fScene("room-and-chair"),
  fRoomX(2.5*m),  // Half-lengths...
  fRoomY(2.5*m),
  fRoomZ(1.3*m),
  fWindowX(10*cm),
  fWindowY(75*cm),
  fWindowZ(50*cm),
  fWindowSillHeight(80*cm),
  fWindowOffset(-50*cm),
  fDoorFrameX(10*cm),
  fDoorFrameY(50*cm),
  fDoorFrameZ(1*m),
  fDoorFrameOffset(1.5*m),
  fDoorX(2*cm),
  fDoorY(50*cm),
  fDoorZ(1*m),
  fChairX(20*cm),          // Half overall width.
  fChairY(20*cm),          // Half overall depth.
  fChairZ(45*cm),          // Half overall height.
  fChairSeat(20*cm),       // Half height of top of seat.
  fChairThickness(3.*cm)   // Half thicknes of back, seat, legs.
{
  new PerspectiveVisActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PerspectiveVisAction::Draw()
{
  fpVisManager = G4VVisManager::GetConcreteInstance();
  if (fpVisManager) {

    // All scenes assume upvector z and origin on "floor"...

    if (fScene == "room-and-chair" )
      {
        RoomAndChair();
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PerspectiveVisAction::RoomAndChair()
{
  // Simple box, size of a room, translated so origin is on xy "floor"...
  G4VisAttributes room_visAtts(G4Colour::Red());
  room_visAtts.SetForceWireframe(true);
  ExtendedDraw
    (G4Box("box",fRoomX,fRoomY,fRoomZ), room_visAtts, G4TranslateZ3D(fRoomZ));

  // Windows...
  G4Box ("window",fWindowX,fWindowY,fWindowZ);
  ExtendedDraw
    (G4Box("window-x",fWindowX,fWindowY,fWindowZ),
     room_visAtts,
     G4Translate3D(-fRoomX - fWindowX, fWindowOffset, fWindowY + fWindowSillHeight));
  ExtendedDraw
    (G4Box("window-y1",fWindowX,fWindowY,fWindowZ),
     room_visAtts,
     G4Translate3D(0., -fRoomY - fWindowX, fWindowY + fWindowSillHeight) *
     G4RotateZ3D(90.*deg));
  ExtendedDraw
    (G4Box("window-y2",fWindowX,fWindowY,fWindowZ),
     room_visAtts,
     G4Translate3D(0., fRoomY + fWindowX, fWindowY + fWindowSillHeight) *
     G4RotateZ3D(-90.*deg));

  // Door...
  ExtendedDraw
    (G4Box("door-frame",fDoorFrameX,fDoorFrameY,fDoorFrameZ),
     room_visAtts,
     G4Translate3D(-fRoomX - fDoorFrameX, fDoorFrameOffset, fDoorFrameZ));
  ExtendedDraw
    (G4Box("door",fDoorX,fDoorY,fDoorZ),
     room_visAtts,
     G4Translate3D(-fRoomX - fDoorX, fDoorFrameOffset, fDoorZ) *
     G4TranslateY3D(fDoorY) *
     G4RotateZ3D(60.*deg) *
     G4TranslateY3D(-fDoorY));  // Last transform operates first.

  // Chair...
  G4VisAttributes chair_visAtts(G4Colour::Cyan());
  G4Transform3D A = G4RotateZ3D(90.*deg);           // Turn through 90 deg.
  G4Transform3D B = G4RotateY3D(90.*deg);           // Lie down.
  G4Transform3D C = G4RotateZ3D(-20.*deg) ;         // Rotate a little.
  G4Transform3D D = G4TranslateZ3D(fChairY);        // Place on floor.
  G4Transform3D E = G4TranslateY3D(-0.5 * fRoomY);  // Move over to the left...
  G4Transform3D chair_transform = E*D*C*B*A;
  Chair(chair_visAtts, chair_transform);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PerspectiveVisAction::Chair
(const G4VisAttributes& visAtts,
 const G4Transform3D& transform)
{
  // Origin is on floor, z = 0, and in the xy centre...
  ExtendedDraw
    (G4Box("chair-back",fChairX, fChairThickness, fChairZ - fChairSeat),
     visAtts, transform *
     G4Translate3D(0.,-fChairY + fChairThickness, fChairZ + fChairSeat));
  ExtendedDraw
    (G4Box("chair-seat",fChairX, fChairY, fChairThickness),
     visAtts, transform * G4TranslateZ3D(-fChairThickness + 2.* fChairSeat));
  for (int i = -1; i < 2; i+=2) {
    for (int j = -1; j < 2; j+=2) {
      ExtendedDraw
        (G4Box("chair-leg",fChairThickness,
               fChairThickness,
               fChairSeat - fChairThickness),
         visAtts, transform *
         G4Translate3D(i * (fChairX - fChairThickness),
                       j * (fChairY - fChairThickness),
                       fChairSeat - fChairThickness));
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PerspectiveVisAction::ExtendedDraw
(const G4VSolid& solid,
 const G4VisAttributes& visAtts,
 const G4Transform3D& transform)
{
  static const G4double extender = 100.*m;
  static const G4Vector3D x(1,0,0);
  static const G4Vector3D y(0,1,0);
  static const G4Vector3D z(0,0,1);

  // Draw extended edges as requested...
  G4bool any = false, A = false, X = false, Y = false, Z = false;
  if (fOptionString.contains("a")) {A = true; any = true;}
  if (fOptionString.contains("x")) {X = true; any = true;}
  if (fOptionString.contains("y")) {Y = true; any = true;}
  if (fOptionString.contains("z")) {Z = true; any = true;}
  if (any)
    {
      G4Polyhedron* polyhedron = solid.GetPolyhedron();
      G4bool isAuxEdgeVisible = false;  // How do I pick this up???   Can't.
      G4bool notLastFace;
      do {
        G4int n;
        G4Point3D nodes[4];
        notLastFace = polyhedron->GetNextFacet(n, nodes);
        G4bool notLastEdge;
        do {
          G4Point3D v1, v2;
          G4int edgeFlag;
          notLastEdge = polyhedron->GetNextEdge(v1, v2, edgeFlag);
          if (isAuxEdgeVisible || edgeFlag > 0) {
            G4Vector3D v21 = v2 - v1;
            // Check for components of actual edge...
            G4Vector3D v21a = v21;
            v21a.transform(transform);
            // G4cout << "v21a: " << v21a << G4endl;
            using namespace std;
            if (A ||
                (Z && abs(v21a.z()) >
                 sqrt(v21a.x()*v21a.x()+v21a.y()*v21a.y())) ||
                (X && abs(v21a.x()) >
                 sqrt(v21a.y()*v21a.y()+v21a.z()*v21a.z())) ||
                (Y && abs(v21a.y()) >
                 sqrt(v21a.x()*v21a.x()+v21a.x()*v21a.z()))) {
              G4Polyline edge;
              edge.SetVisAttributes(G4Colour(.2,.2,.2));
              edge.push_back(v1 - extender * v21.unit());
              edge.push_back(v2 + extender * v21.unit());
              fpVisManager->Draw(edge, transform);
            }
          }
        } while (notLastEdge);
      } while (notLastFace);
    }

  // Draw actual object...
  fpVisManager->Draw(solid, visAtts, transform);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
