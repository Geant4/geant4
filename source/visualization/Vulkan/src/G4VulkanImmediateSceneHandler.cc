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
// Alex Vergara Gil  17th October 2026
// Ported to Vulkan from the OpenGL work by Andrew Walkden (10th February 1997).
// Vulkan immediate scene - draws immediately to buffer
//                           (saving space on server).

#include "G4VulkanImmediateSceneHandler.hh"

#include "G4AttHolder.hh"
#include "G4Circle.hh"
#include "G4VulkanTransform3D.hh"
#include "G4VulkanViewer.hh"
#include "G4Polyhedron.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Square.hh"
#include "G4Text.hh"

#include <typeinfo>

G4VulkanImmediateSceneHandler::G4VulkanImmediateSceneHandler(G4VGraphicsSystem& system,
                                                             const G4String& name)
  : G4VulkanSceneHandler(system, fSceneIdCount++, name)
{}

G4VulkanImmediateSceneHandler::~G4VulkanImmediateSceneHandler() {}

#include <iomanip>

G4bool G4VulkanImmediateSceneHandler::AddPrimitivePreamble(const G4VMarker& visible)
{
  return AddPrimitivePreambleInternal(visible, true, false);
}
G4bool G4VulkanImmediateSceneHandler::AddPrimitivePreamble(const G4Polyline& visible)
{
  return AddPrimitivePreambleInternal(visible, false, true);
}
G4bool G4VulkanImmediateSceneHandler::AddPrimitivePreamble(const G4Polyhedron& visible)
{
  return AddPrimitivePreambleInternal(visible, false, false);
}

G4bool G4VulkanImmediateSceneHandler::AddPrimitivePreambleInternal(const G4Visible& visible,
                                                                   bool isMarker, bool isPolyline)
{
  // Get applicable vis attributes for all primitives.
  fpVisAttribs = fpViewer->GetApplicableVisAttributes(visible.GetVisAttributes());
  const G4Colour& c = GetColour();
  G4double opacity = c.GetAlpha();

  G4bool transparency_enabled = true;
  G4bool isMarkerNotHidden = true;
  G4VulkanViewer* pViewer = dynamic_cast<G4VulkanViewer*>(fpViewer);
  if (pViewer)
  {
    transparency_enabled = pViewer->transparency_enabled;
    isMarkerNotHidden = pViewer->fVP.IsMarkerNotHidden();
  }

  G4bool isMarkerOrPolyline = isMarker || isPolyline;
  G4bool treatAsTransparent = transparency_enabled && opacity < 1.;
  G4bool treatAsNotHidden = isMarkerNotHidden && (isMarker || isPolyline);

  if (fProcessing2D)
    vkDisable(VK_DEPTH_TEST);
  else
  {
    if (isMarkerOrPolyline && isMarkerNotHidden)
      vkDisable(VK_DEPTH_TEST);
    else
    {
      vkEnable(VK_DEPTH_TEST);
      vkDepthFunc(VK_LEQUAL);
    }
  }

  if (fThreePassCapable)
  {
    // Ensure transparent objects are drawn opaque ones and before
    // non-hidden markers.  The problem of blending/transparency/alpha
    // is quite a tricky one - see History of opengl-V07-01-01/2/3.
    if (!(fSecondPassForTransparency || fThirdPassForNonHiddenMarkers))
    {
      // First pass...
      if (treatAsTransparent)
      {  // Request pass for transparent objects...
        fSecondPassForTransparencyRequested = true;
      }
      if (treatAsNotHidden)
      {  // Request pass for non-hidden markers...
        fThirdPassForNonHiddenMarkersRequested = true;
      }
      // On first pass, transparent objects and non-hidden markers are not drawn...
      if (treatAsTransparent || treatAsNotHidden)
      {
        return false;
      }
    }

    // On second pass, only transparent objects are drawn...
    if (fSecondPassForTransparency)
    {
      if (!treatAsTransparent)
      {
        return false;
      }
    }

    // On third pass, only non-hidden markers are drawn...
    if (fThirdPassForNonHiddenMarkers)
    {
      if (!treatAsNotHidden)
      {
        return false;
      }
    }
  }  // fThreePassCapable

  // Loads G4Atts for picking...
  if (fpViewer->GetViewParameters().IsPicking())
  {
    vkLoadName(++fPickName);
    G4AttHolder* holder = new G4AttHolder;
    LoadAtts(visible, holder);
    fPickMap[fPickName] = holder;
  }

  if (transparency_enabled)
  {
    vkColor4d(c.GetRed(), c.GetGreen(), c.GetBlue(), c.GetAlpha());
  }
  else
  {
    vkColor3d(c.GetRed(), c.GetGreen(), c.GetBlue());
  }

  return true;
}

void G4VulkanImmediateSceneHandler::AddPrimitive(const G4Polyline& polyline)
{
  G4bool furtherprocessing = AddPrimitivePreamble(polyline);
  if (furtherprocessing)
  {
    G4VulkanSceneHandler::AddPrimitive(polyline);
  }
}

void G4VulkanImmediateSceneHandler::AddPrimitive(const G4Polymarker& polymarker)
{
  G4bool furtherprocessing = AddPrimitivePreamble(polymarker);
  if (furtherprocessing)
  {
    G4VulkanSceneHandler::AddPrimitive(polymarker);
  }
}

void G4VulkanImmediateSceneHandler::AddPrimitive(const G4Text& text)
{
  // Note: colour is still handled in
  // G4VulkanSceneHandler::AddPrimitive(const G4Text&).
  G4bool furtherprocessing = AddPrimitivePreamble(text);
  if (furtherprocessing)
  {
    G4VulkanSceneHandler::AddPrimitive(text);
  }
}

void G4VulkanImmediateSceneHandler::AddPrimitive(const G4Circle& circle)
{
  G4bool furtherprocessing = AddPrimitivePreamble(circle);
  if (furtherprocessing)
  {
    G4VulkanSceneHandler::AddPrimitive(circle);
  }
}

void G4VulkanImmediateSceneHandler::AddPrimitive(const G4Square& square)
{
  G4bool furtherprocessing = AddPrimitivePreamble(square);
  if (furtherprocessing)
  {
    G4VulkanSceneHandler::AddPrimitive(square);
  }
}

void G4VulkanImmediateSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron)
{
  // Note: colour is still handled in
  // G4VulkanSceneHandler::AddPrimitive(const G4Polyhedron&).
  G4bool furtherprocessing = AddPrimitivePreamble(polyhedron);
  if (furtherprocessing)
  {
    G4VulkanSceneHandler::AddPrimitive(polyhedron);
  }
}

void G4VulkanImmediateSceneHandler::BeginPrimitives(const G4Transform3D& objectTransformation)
{
  G4VulkanSceneHandler::BeginPrimitives(objectTransformation);

  G4VulkanTransform3D vkt(objectTransformation);

  vkPushMatrix();

  /*************************** Check matrix.
  const VKdouble* m = vkt.GetVKMatrix ();
  G4cout << "G4VulkanTransform3D matrix:";
  for (int i = 0; i < 16; i++) {
    if ((i % 4) == 0) G4cout << '\n';
    G4cout << std::setw (15) << m[i];
  }
  G4cout << G4endl;
  *****************************************/

  vkMultMatrixd(vkt.GetVKMatrix());
}

void G4VulkanImmediateSceneHandler::EndPrimitives()
{
  vkPopMatrix();

  // See all primitives immediately...  At least soon...
  ScaledFlush();

  G4VulkanSceneHandler::EndPrimitives();
}

void G4VulkanImmediateSceneHandler::BeginPrimitives2D(const G4Transform3D& objectTransformation)
{
  G4VulkanSceneHandler::BeginPrimitives2D(objectTransformation);

  // Push current 3D world matrices and load identity to define screen
  // coordinates...
  vkMatrixMode(VK_PROJECTION);
  vkPushMatrix();
  vkLoadIdentity();
  G4VulkanViewer* pViewer = dynamic_cast<G4VulkanViewer*>(fpViewer);
  if (pViewer)
  {
    pViewer->g4VkOrtho(-1., 1., -1., 1., -G4VULKAN_FLT_BIG, G4VULKAN_FLT_BIG);
  }
  vkMatrixMode(VK_MODELVIEW);
  vkPushMatrix();
  vkLoadIdentity();
  G4VulkanTransform3D vkt(objectTransformation);
  vkMultMatrixd(vkt.GetVKMatrix());
  vkDisable(VK_DEPTH_TEST);  // But see parent scene handler!!  In
  vkDisable(VK_LIGHTING);  // some cases, we need to re-iterate this.
}

void G4VulkanImmediateSceneHandler::EndPrimitives2D()
{
  // Pop current 3D world matrices back again...
  vkMatrixMode(VK_PROJECTION);
  vkPopMatrix();
  vkMatrixMode(VK_MODELVIEW);
  vkPopMatrix();

  // See all primitives immediately...  At least soon...
  ScaledFlush();

  G4VulkanSceneHandler::EndPrimitives2D();
}

void G4VulkanImmediateSceneHandler::BeginModeling()
{
  G4VSceneHandler::BeginModeling();
}

void G4VulkanImmediateSceneHandler::EndModeling()
{
  G4VSceneHandler::EndModeling();
}

void G4VulkanImmediateSceneHandler::ClearTransientStore()
{
  // Nothing to do except redraw the scene ready for the next event.
  if (fpViewer)
  {
    fpViewer->SetView();
    fpViewer->ClearView();
    fpViewer->DrawView();
  }
}

G4int G4VulkanImmediateSceneHandler::fSceneIdCount = 0;
