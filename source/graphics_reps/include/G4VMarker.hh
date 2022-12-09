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
// G4VMarker - base class for markers - circles, squares, etc.
// John Allison  17/11/96.

// Class Description:
// G4VMarkers are 2-dimensional G4VVisPrims with the special
// properties (a) of always facing the camera and (b) of having the
// possibility of a size defined in screen units (pixels) or paper
// units (points - there are 72 points per inch).  The convention is
// that if a world size is not specified, then the marker will be
// drawn to the given screen size or paper size independent of the
// viewing transformation in effect.
//
// "Size" means "overall size", e.g., diameter of circle, side of
// square, height of text (but diameter and radius access functions
// are defined to avoid ambiguity).
//
// So the user who constructs the marker decides whether it should be
// drawn to a given size in world coordinates by setting the world
// size.  Alternatively, the user can set the screen size (internally,
// the world size is set to zero) and the marker is drawn to its
// screen size.  Finally, the user may decide not to set any size; in
// that case, it is drawn according to the sizes specified in the
// default marker specified in G4ViewParameters.
//
// Also in G4ViewParameters is a "global marker scale" which is a
// factor by which all marker sizes are multiplied before drawing.
//
// Thus the graphics system driver scene handler code might look like:
//
// void G4XXXGraphicsSceneHandler::AddPrimitive (const G4Circle& circle) {
//   G4bool hidden = !(fpView -> GetViewParameters().IsMarkerNotHidden());
//   const G4Colour&     colour = GetColour (circle);  
//                                   // Base class GetColour.
//   G4VMarker::FillStyle style = circle.GetFillStyle();
//   const G4Point3D&    centre = circle.GetPosition();
//   MarkerSizeType sizeType;
//   G4double size = GetMarkerSize (circle, sizeType);
//   switch (sizeType) {
//   default:
//   case screen:
//     // Draw in screen coordinates.
//     // ...
//     break;
//   case world:
//     // Draw in world coordinates.
//     // ...
//     break;
//   }
// }
// Class Description - End:


#ifndef G4VMARKER_HH
#define G4VMARKER_HH

#include "globals.hh"
#include "G4Visible.hh"
#include "G4Point3D.hh"
#include "G4Colour.hh"
#include "G4Color.hh"

class G4VMarker: public G4Visible {

  friend std::ostream& operator << (std::ostream& os, const G4VMarker&);

public: // With description

  enum FillStyle {noFill, hashed, filled};
  enum SizeType {none, world, screen};

  //////////////////////////////////////////////////////
  // Constructors...
  G4VMarker ();
  G4VMarker (const G4VMarker&) = default;
  G4VMarker (G4VMarker&& ) = default;
  G4VMarker (const G4Point3D& position);

  //////////////////////////////////////////////////////
  // Destructor...
  ~G4VMarker () override;

  //////////////////////////////////////////////////////
  // Assignment...
  G4VMarker& operator = (const G4VMarker&) = default;
  G4VMarker& operator = (G4VMarker&&) = default;

  //////////////////////////////////////////////////////
  // Logical...
  G4bool operator != (const G4VMarker&) const;
  G4bool operator == (const G4VMarker&) const;

  /////////////////////////////////////////////////////
  // Get functions...
  G4Point3D GetPosition       () const;
  SizeType  GetSizeType       () const;
  G4double  GetWorldSize      () const;
  G4double  GetWorldDiameter  () const;
  G4double  GetWorldRadius    () const;
  G4double  GetScreenSize     () const;
  G4double  GetScreenDiameter () const;
  G4double  GetScreenRadius   () const;
  FillStyle GetFillStyle      () const;

  /////////////////////////////////////////////////////
  // Set functions...
  void SetPosition       (const G4Point3D&);
  void SetSize           (SizeType, G4double);
  void SetDiameter       (SizeType, G4double);
  void SetRadius         (SizeType, G4double);
  void SetWorldSize      (G4double);
  void SetWorldDiameter  (G4double);
  void SetWorldRadius    (G4double);
  void SetScreenSize     (G4double);
  void SetScreenDiameter (G4double);
  void SetScreenRadius   (G4double);
  void SetFillStyle      (FillStyle);

private:
  G4Point3D fPosition;
  G4double  fWorldSize;   // Default 0. means use screen size.
  G4double  fScreenSize;  // Default 0. means use global default.
  FillStyle fFillStyle;
};

#include "G4VMarker.icc"

#endif
