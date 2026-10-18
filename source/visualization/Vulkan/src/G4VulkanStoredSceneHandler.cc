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
// Vulkan stored scene - creates Vulkan display lists.

#include "G4VulkanStoredSceneHandler.hh"

#include "G4AttHolder.hh"
#include "G4Circle.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeModel.hh"
#include "G4VulkanTransform3D.hh"
#include "G4VulkanViewer.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4Polyhedron.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Square.hh"
#include "G4Text.hh"
#include "G4VPhysicalVolume.hh"

#include <typeinfo>

G4int G4VulkanStoredSceneHandler::fSceneIdCount = 0;

G4int G4VulkanStoredSceneHandler::fDisplayListId = 0;

G4VulkanStoredSceneHandler::PO::PO()
  : fDisplayListId(0), fPickName(0), fpG4TextPlus(0), fMarkerOrPolyline(false)
{}

G4VulkanStoredSceneHandler::PO::PO(const G4VulkanStoredSceneHandler::PO& po)
  : fDisplayListId(po.fDisplayListId),
    fTransform(po.fTransform),
    fPickName(po.fPickName),
    fColour(po.fColour),
    fpG4TextPlus(po.fpG4TextPlus ? new G4TextPlus(*po.fpG4TextPlus) : 0),
    fMarkerOrPolyline(po.fMarkerOrPolyline)
{}

G4VulkanStoredSceneHandler::PO::PO(G4int id, const G4Transform3D& tr)
  : fDisplayListId(id), fTransform(tr), fPickName(0), fpG4TextPlus(0), fMarkerOrPolyline(false)
{}

G4VulkanStoredSceneHandler::PO::~PO()
{
  delete fpG4TextPlus;
}

G4VulkanStoredSceneHandler::PO&
G4VulkanStoredSceneHandler::PO::operator=(const G4VulkanStoredSceneHandler::PO& rhs)
{
  if (&rhs == this) return *this;
  fDisplayListId = rhs.fDisplayListId;
  fTransform = rhs.fTransform;
  fPickName = rhs.fPickName;
  fColour = rhs.fColour;
  fpG4TextPlus = rhs.fpG4TextPlus ? new G4TextPlus(*rhs.fpG4TextPlus) : 0;
  fMarkerOrPolyline = rhs.fMarkerOrPolyline;
  return *this;
}

G4VulkanStoredSceneHandler::TO::TO()
  : fDisplayListId(0),
    fPickName(0),
    fStartTime(-G4VisAttributes::fVeryLongTime),
    fEndTime(G4VisAttributes::fVeryLongTime),
    fpG4TextPlus(0),
    fMarkerOrPolyline(false)
{}

G4VulkanStoredSceneHandler::TO::TO(const G4VulkanStoredSceneHandler::TO& to)
  : fDisplayListId(to.fDisplayListId),
    fTransform(to.fTransform),
    fPickName(to.fPickName),
    fStartTime(to.fStartTime),
    fEndTime(to.fEndTime),
    fColour(to.fColour),
    fpG4TextPlus(to.fpG4TextPlus ? new G4TextPlus(*to.fpG4TextPlus) : 0),
    fMarkerOrPolyline(to.fMarkerOrPolyline)
{}

G4VulkanStoredSceneHandler::TO::TO(G4int id, const G4Transform3D& tr)
  : fDisplayListId(id),
    fTransform(tr),
    fPickName(0),
    fStartTime(-G4VisAttributes::fVeryLongTime),
    fEndTime(G4VisAttributes::fVeryLongTime),
    fpG4TextPlus(0),
    fMarkerOrPolyline(false)
{}

G4VulkanStoredSceneHandler::TO::~TO()
{
  delete fpG4TextPlus;
}

G4VulkanStoredSceneHandler::TO&
G4VulkanStoredSceneHandler::TO::operator=(const G4VulkanStoredSceneHandler::TO& rhs)
{
  if (&rhs == this) return *this;
  fDisplayListId = rhs.fDisplayListId;
  fTransform = rhs.fTransform;
  fPickName = rhs.fPickName;
  fStartTime = rhs.fStartTime;
  fEndTime = rhs.fEndTime;
  fColour = rhs.fColour;
  fpG4TextPlus = rhs.fpG4TextPlus ? new G4TextPlus(*rhs.fpG4TextPlus) : 0;
  fMarkerOrPolyline = rhs.fMarkerOrPolyline;
  return *this;
}

G4VulkanStoredSceneHandler::G4VulkanStoredSceneHandler(G4VGraphicsSystem& system,
                                                       const G4String& name)
  : G4VulkanSceneHandler(system, fSceneIdCount++, name), fDoNotUseDisplayList(false), fTopPODL(0)
{}

G4VulkanStoredSceneHandler::~G4VulkanStoredSceneHandler() {}

void G4VulkanStoredSceneHandler::BeginPrimitives(const G4Transform3D& objectTransformation)
{
  G4VulkanSceneHandler::BeginPrimitives(objectTransformation);
  if (fReadyForTransients) vkDrawBuffer(VK_FRONT);
  // Display list setup moved to AddPrimitivePreamble.  See notes there.
}

void G4VulkanStoredSceneHandler::EndPrimitives()
{
  // See all primitives immediately...  At least soon...
  ScaledFlush();
  vkDrawBuffer(VK_BACK);
  G4VulkanSceneHandler::EndPrimitives();
}

void G4VulkanStoredSceneHandler::BeginPrimitives2D(const G4Transform3D& objectTransformation)
{
  G4VulkanSceneHandler::BeginPrimitives2D(objectTransformation);
  if (fReadyForTransients) vkDrawBuffer(VK_FRONT);
}

void G4VulkanStoredSceneHandler::EndPrimitives2D()
{
  // See all primitives immediately...  At least soon...
  ScaledFlush();
  vkDrawBuffer(VK_BACK);
  G4VulkanSceneHandler::EndPrimitives2D();
}

G4bool G4VulkanStoredSceneHandler::AddPrimitivePreamble(const G4VMarker& visible)
{
  return AddPrimitivePreambleInternal(visible, true, false);
}
G4bool G4VulkanStoredSceneHandler::AddPrimitivePreamble(const G4Polyline& visible)
{
  return AddPrimitivePreambleInternal(visible, false, true);
}
G4bool G4VulkanStoredSceneHandler::AddPrimitivePreamble(const G4Polyhedron& visible)
{
  return AddPrimitivePreambleInternal(visible, false, false);
}

G4bool G4VulkanStoredSceneHandler::AddPrimitivePreambleInternal(const G4Visible& visible,
                                                                bool isMarker, bool isPolyline)
{
  // Get applicable vis attributes for all primitives.
  fpVisAttribs = fpViewer->GetApplicableVisAttributes(visible.GetVisAttributes());
  const G4Colour& c = GetColour();
  G4double opacity = c.GetAlpha();

  G4bool transparency_enabled = true;
  G4bool isMarkerNotHidden = true;
  G4VulkanViewer* pVKViewer = dynamic_cast<G4VulkanViewer*>(fpViewer);
  if (pVKViewer)
  {
    transparency_enabled = pVKViewer->transparency_enabled;
    isMarkerNotHidden = pVKViewer->fVP.IsMarkerNotHidden();
  }

  G4bool isTransparent = opacity < 1.;
  G4bool isMarkerOrPolyline = isMarker || isPolyline;
  G4bool treatAsTransparent = transparency_enabled && isTransparent;
  G4bool treatAsNotHidden = isMarkerNotHidden && isMarkerOrPolyline;

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
    // Ensure transparent objects are drawn *after* opaque ones and before
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
        return false;  // No further processing.
      }
    }

    // On second pass, only transparent objects are drawn...
    if (fSecondPassForTransparency)
    {
      if (!treatAsTransparent)
      {
        return false;  // No further processing.
      }
    }

    // On third pass, only non-hidden markers are drawn...
    if (fThirdPassForNonHiddenMarkers)
    {
      if (!treatAsNotHidden)
      {
        return false;  // No further processing.
      }
    }
  }  // fThreePassCapable

  // Loads G4Atts for picking...
  G4bool isPicking = false;
  if (fpViewer->GetViewParameters().IsPicking())
  {
    isPicking = true;
    vkLoadName(++fPickName);
    G4AttHolder* holder = new G4AttHolder;
    LoadAtts(visible, holder);
    fPickMap[fPickName] = holder;
  }

  // Because of our need to control colour of transients (display by
  // time fading), display lists may only cover a single primitive.
  // So display list setup is here.

  if (fDoNotUseDisplayList)
  {
    vkPushMatrix();
    G4VulkanTransform3D vkt(fObjectTransformation);
    vkMultMatrixd(vkt.GetVKMatrix());
    if (transparency_enabled)
    {
      vkColor4d(c.GetRed(), c.GetGreen(), c.GetBlue(), c.GetAlpha());
    }
    else
    {
      vkColor3d(c.GetRed(), c.GetGreen(), c.GetBlue());
    }
  }
  else
  {
    fDisplayListId = vkGenLists(1);
    if (vkGetError() == VK_OUT_OF_MEMORY)
    {
      static G4int errorCount = 0;
      if (errorCount < 5)
      {
        errorCount++;
        G4ExceptionDescription ed;
        ed << "Error attempting to create an Vulkan display list."
              "\nCurrent display list id: "
           << fDisplayListId << "\nMaybe out of memory?";
        G4Exception("G4VulkanStoredSceneHandler::AddPrimitivePreambleInternal", "vulkan1001",
                    JustWarning, ed);
      }
      return false;
    }
    if (fReadyForTransients)
    {
      TO to(fDisplayListId, fObjectTransformation);
      if (isPicking) to.fPickName = fPickName;
      to.fColour = c;
      to.fStartTime = fpVisAttribs->GetStartTime();
      to.fEndTime = fpVisAttribs->GetEndTime();
      to.fMarkerOrPolyline = isMarkerOrPolyline;
      fTOList.push_back(to);
      // For transient objects, colour, transformation, are kept in
      // the TO, so should *not* be in the display list.  As mentioned
      // above, in some cases (display-by-time fading) we need to have
      // independent control of colour.  But for now transform and set
      // colour for immediate display.
      vkPushMatrix();
      G4VulkanTransform3D vkt(fObjectTransformation);
      vkMultMatrixd(vkt.GetVKMatrix());
      if (transparency_enabled)
      {
        vkColor4d(c.GetRed(), c.GetGreen(), c.GetBlue(), c.GetAlpha());
      }
      else
      {
        vkColor3d(c.GetRed(), c.GetGreen(), c.GetBlue());
      }
      (void)ExtraTOProcessing(visible, fTOList.size() - 1);
      // Ignore return value of the above.  If this visible does not use
      // gl commands, a display list is created that is empty and not
      // used.
      vkNewList(fDisplayListId, VK_COMPILE_AND_EXECUTE);
    }
    else
    {
      PO po(fDisplayListId, fObjectTransformation);
      if (isPicking) po.fPickName = fPickName;
      po.fColour = c;
      po.fMarkerOrPolyline = isMarkerOrPolyline;
      fPOList.push_back(po);
      // For permanent objects, colour is kept in the PO, so should
      // *not* be in the display list.  This is so that sub-classes
      // may implement colour modifications according to their own
      // criteria, e.g., scene tree slider in Qt.  But for now set
      // colour for immediate display.
      if (transparency_enabled)
      {
        vkColor4d(c.GetRed(), c.GetGreen(), c.GetBlue(), c.GetAlpha());
      }
      else
      {
        vkColor3d(c.GetRed(), c.GetGreen(), c.GetBlue());
      }
      G4bool usesVKCommands = ExtraPOProcessing(visible, fPOList.size() - 1);
      // Transients are displayed as they come (VK_COMPILE_AND_EXECUTE
      // above) but persistents are compiled into display lists
      // (VK_COMPILE only) and then drawn from the display lists with
      // their fObjectTransformation as stored in fPOList.  Thus,
      // there is no need to do vkMultMatrixd here.  If
      // ExtraPOProcessing says the visible object does not use gl
      // commands, simply return and abandon further processing.  It
      // is assumed that all relevant information is kept in the
      // POList.
      if (!usesVKCommands) return false;
      vkNewList(fDisplayListId, VK_COMPILE);
    }
  }

  if (fProcessing2D)
  {
    // Push current 3D world matrices and load identity to define screen
    // coordinates...
    vkMatrixMode(VK_PROJECTION);
    vkPushMatrix();
    vkLoadIdentity();
    if (pVKViewer)
    {
      pVKViewer->g4VkOrtho(-1., 1., -1., 1., -G4VULKAN_FLT_BIG, G4VULKAN_FLT_BIG);
    }
    vkMatrixMode(VK_MODELVIEW);
    vkPushMatrix();
    vkLoadIdentity();
    G4VulkanTransform3D vkt(fObjectTransformation);
    vkMultMatrixd(vkt.GetVKMatrix());
    vkDisable(VK_LIGHTING);
  }
  else
  {
    if (isMarker)
    {
      vkDisable(VK_LIGHTING);
    }
    else
    {
      vkEnable(VK_LIGHTING);
    }
  }

  return true;
}

void G4VulkanStoredSceneHandler::AddPrimitivePostamble()
{
  if (fProcessing2D)
  {
    // Pop current 3D world matrices back again...
    vkMatrixMode(VK_PROJECTION);
    vkPopMatrix();
    vkMatrixMode(VK_MODELVIEW);
    vkPopMatrix();
  }

  //  if ((vkGetError() == VK_TABLE_TOO_LARGE) || (vkGetError() == VK_OUT_OF_MEMORY)) {  // Could
  //  close?
  if (vkGetError() == VK_OUT_OF_MEMORY)
  {  // Could close?
    G4cerr << "ERROR: G4VulkanStoredSceneHandler::AddPrimitivePostamble: Failure"
              "  to allocate display List for fTopPODL - try Vulkan Immediated mode."
           << G4endl;
  }
  if (!fDoNotUseDisplayList)
  {
    vkEndList();
    if (vkGetError() == VK_OUT_OF_MEMORY)
    {  // Could close?
      G4cerr << "ERROR: G4VulkanStoredSceneHandler::AddPrimitivePostamble: Failure"
                "  to allocate display List for fTopPODL - try Vulkan Immediated mode."
             << G4endl;
    }
  }
  if (fReadyForTransients || fDoNotUseDisplayList)
  {
    vkPopMatrix();
  }
}

void G4VulkanStoredSceneHandler::AddPrimitive(const G4Polyline& polyline)
{
  G4bool furtherprocessing = AddPrimitivePreamble(polyline);
  if (furtherprocessing)
  {
    G4VulkanSceneHandler::AddPrimitive(polyline);
    AddPrimitivePostamble();
  }
}

void G4VulkanStoredSceneHandler::AddPrimitive(const G4Polymarker& polymarker)
{
  G4bool furtherprocessing = AddPrimitivePreamble(polymarker);
  if (furtherprocessing)
  {
    G4VulkanSceneHandler::AddPrimitive(polymarker);
    AddPrimitivePostamble();
  }
}

void G4VulkanStoredSceneHandler::AddPrimitive(const G4Text& text)
{
  // Note: colour is still handled in
  // G4VulkanSceneHandler::AddPrimitive(const G4Text&), so it still
  // gets into the display list
  G4bool furtherprocessing = AddPrimitivePreamble(text);
  if (furtherprocessing)
  {
    G4VulkanSceneHandler::AddPrimitive(text);
    AddPrimitivePostamble();
  }
}

void G4VulkanStoredSceneHandler::AddPrimitive(const G4Circle& circle)
{
  G4bool furtherprocessing = AddPrimitivePreamble(circle);
  if (furtherprocessing)
  {
    G4VulkanSceneHandler::AddPrimitive(circle);
    AddPrimitivePostamble();
  }
}

void G4VulkanStoredSceneHandler::AddPrimitive(const G4Square& square)
{
  G4bool furtherprocessing = AddPrimitivePreamble(square);
  if (furtherprocessing)
  {
    G4VulkanSceneHandler::AddPrimitive(square);
    AddPrimitivePostamble();
  }
}

void G4VulkanStoredSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron)
{
  // Note: colour is still handled in
  // G4VulkanSceneHandler::AddPrimitive(const G4Polyhedron&), so it still
  // gets into the display list
  G4bool furtherprocessing = AddPrimitivePreamble(polyhedron);
  if (furtherprocessing)
  {
    G4VulkanSceneHandler::AddPrimitive(polyhedron);
    AddPrimitivePostamble();
  }
}

void G4VulkanStoredSceneHandler::BeginModeling()
{
  G4VSceneHandler::BeginModeling();
  /* Debug...
  fDisplayListId = vkGenLists (1);
  G4cout << "VK::fDisplayListId (start): " << fDisplayListId << G4endl;
  */
}

void G4VulkanStoredSceneHandler::EndModeling()
{
  // Make a List which calls the other lists.
  fTopPODL = vkGenLists(1);
  if (vkGetError() == VK_OUT_OF_MEMORY)
  {  // Could pre-allocate?
    G4cerr << "ERROR: G4VulkanStoredSceneHandler::EndModeling: Failure to allocate"
              "  display List for fTopPODL - try Vulkan Immediated mode."
           << G4endl;
  }
  else
  {
    vkNewList(fTopPODL, VK_COMPILE);
    {
      for (size_t i = 0; i < fPOList.size(); i++)
      {
        vkPushMatrix();
        G4VulkanTransform3D vkt(fPOList[i].fTransform);
        vkMultMatrixd(vkt.GetVKMatrix());
        if (fpViewer->GetViewParameters().IsPicking()) vkLoadName(fPOList[i].fPickName);
        vkCallList(fPOList[i].fDisplayListId);
        vkPopMatrix();
      }
    }
    vkEndList();

    if (vkGetError() == VK_OUT_OF_MEMORY)
    {  // Could close?
      G4cerr << "ERROR: G4VulkanStoredSceneHandler::EndModeling: Failure to allocate"
                "  display List for fTopPODL - try Vulkan Immediated mode."
             << G4endl;
    }
  }

  G4VSceneHandler::EndModeling();
}

void G4VulkanStoredSceneHandler::ClearStore()
{
  // G4cout << "G4VulkanStoredSceneHandler::ClearStore" << G4endl;

  G4VSceneHandler::ClearStore();  // Sets need kernel visit, etc.

  // Delete Vulkan permanent display lists.
  for (size_t i = 0; i < fPOList.size(); i++)
    vkDeleteLists(fPOList[i].fDisplayListId, 1);
  if (fTopPODL) vkDeleteLists(fTopPODL, 1);
  fTopPODL = 0;

  // Clear other lists, dictionary, etc.
  fPOList.clear();
  fSolidMap.clear();
  ClearAndDestroyAtts();

  // ...and clear transient store...
  for (size_t i = 0; i < fTOList.size(); i++)
    vkDeleteLists(fTOList[i].fDisplayListId, 1);
  fTOList.clear();
}

void G4VulkanStoredSceneHandler::ClearTransientStore()
{
  // G4cout << "G4VulkanStoredSceneHandler::ClearTransientStore" << G4endl;

  // Delete Vulkan transient display lists and Transient Objects themselves.
  for (size_t i = 0; i < fTOList.size(); i++)
    vkDeleteLists(fTOList[i].fDisplayListId, 1);
  fTOList.clear();

  // Redraw the scene ready for the next event.
  if (fpViewer)
  {
    fpViewer->SetView();
    fpViewer->ClearView();
    fpViewer->DrawView();
  }
}
