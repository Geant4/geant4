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
// Ported to Vulkan from the OpenGL work by Andrew Walkden (7th February 1997).
// Class G4VulkanStoredViewer : Encapsulates the `storedness' of
//                            an Vulkan view, for inheritance by
//                            derived (X, Xm...) classes.

#include "G4VulkanStoredViewer.hh"

#include "G4Circle.hh"
#include "G4VulkanStoredSceneHandler.hh"
#include "G4VulkanTransform3D.hh"
#include "G4PhysicalConstants.hh"
#include "G4Scene.hh"
#include "G4Text.hh"
#include "G4UnitsTable.hh"

G4VulkanStoredViewer::G4VulkanStoredViewer(G4VulkanStoredSceneHandler& sceneHandler)
  : G4VViewer(sceneHandler, -1),
    G4VulkanViewer(sceneHandler),
    fG4VulkanStoredSceneHandler(sceneHandler),
    fDepthTestEnable(true)
{
  fLastVP = fDefaultVP;  // Update in sub-class after KernelVisitDecision
}

G4VulkanStoredViewer::~G4VulkanStoredViewer() {}

void G4VulkanStoredViewer::KernelVisitDecision()
{
  // If there's a significant difference with the last view parameters
  // of either the scene handler or this viewer, trigger a rebuild.

  if (!fG4VulkanStoredSceneHandler.fTopPODL || CompareForKernelVisit(fLastVP))
  {
    NeedKernelVisit();
  }
}

G4bool G4VulkanStoredViewer::CompareForKernelVisit(G4ViewParameters& lastVP)
{
  if ((lastVP.GetDrawingStyle() != fVP.GetDrawingStyle())
      || (lastVP.GetNumberOfCloudPoints() != fVP.GetNumberOfCloudPoints())
      || (lastVP.IsAuxEdgeVisible() != fVP.IsAuxEdgeVisible())
      || (lastVP.IsCulling() != fVP.IsCulling())
      || (lastVP.IsCullingInvisible() != fVP.IsCullingInvisible())
      || (lastVP.IsDensityCulling() != fVP.IsDensityCulling())
      || (lastVP.IsCullingCovered() != fVP.IsCullingCovered())
      || (lastVP.GetCBDAlgorithmNumber() != fVP.GetCBDAlgorithmNumber()) ||
      // Note: Section and Cutaway can reveal back-facing faces. If
      // backface culling is implemented, the image can look strange because
      // the back-facing faces are not there. For the moment, we have disabled
      // (commented out) backface culling (it seems not to affect performance -
      // in fact, performance seems to improve), so there is no problem.
      (lastVP.IsSection() != fVP.IsSection()) ||
      // Section (DCUT) is NOT implemented locally so we need to visit the kernel.
      // (lastVP.IsCutaway ()          != fVP.IsCutaway ())          ||
      // Cutaways are implemented locally so we do not need to visit the kernel.
      (lastVP.IsExplode() != fVP.IsExplode()) || (lastVP.GetNoOfSides() != fVP.GetNoOfSides())
      || (lastVP.GetGlobalMarkerScale() != fVP.GetGlobalMarkerScale())
      || (lastVP.GetGlobalLineWidthScale() != fVP.GetGlobalLineWidthScale())
      || (lastVP.IsMarkerNotHidden() != fVP.IsMarkerNotHidden())
      || (lastVP.GetDefaultVisAttributes()->GetColour()
          != fVP.GetDefaultVisAttributes()->GetColour())
      || (lastVP.GetDefaultTextVisAttributes()->GetColour()
          != fVP.GetDefaultTextVisAttributes()->GetColour())
      || (lastVP.GetBackgroundColour() != fVP.GetBackgroundColour())
      || (lastVP.IsPicking() != fVP.IsPicking())
      || (lastVP.GetVisAttributesModifiers() != fVP.GetVisAttributesModifiers())
      || (lastVP.IsSpecialMeshRendering() != fVP.IsSpecialMeshRendering())
      || (lastVP.GetSpecialMeshRenderingOption() != fVP.GetSpecialMeshRenderingOption())
      || (lastVP.GetTransparencyByDepth() != fVP.GetTransparencyByDepth())
      || (lastVP.IsDotsSmooth() != fVP.IsDotsSmooth())
      || (lastVP.GetDotsSize() != fVP.GetDotsSize()))
    return true;

  if (lastVP.IsDensityCulling() && (lastVP.GetVisibleDensity() != fVP.GetVisibleDensity()))
    return true;

  //  /**************************************************************
  //   If section (DCUT) is implemented locally, comment this out.
  if (lastVP.IsSection() && (lastVP.GetSectionPlane() != fVP.GetSectionPlane())) return true;
  //   ***************************************************************/

  /**************************************************************
   If cutaways are implemented locally, comment this out.
   if (lastVP.IsCutaway ()) {
   if (vp.GetCutawayMode() != fVP.GetCutawayMode()) return true;
   if (lastVP.GetCutawayPlanes ().size () !=
   fVP.GetCutawayPlanes ().size ()) return true;
   for (size_t i = 0; i < lastVP.GetCutawayPlanes().size(); ++i)
   if (lastVP.GetCutawayPlanes()[i] != fVP.GetCutawayPlanes()[i])
   return true;
   }
   ***************************************************************/

  if (lastVP.GetCBDAlgorithmNumber() > 0)
  {
    if (lastVP.GetCBDParameters().size() != fVP.GetCBDParameters().size())
      return true;
    else if (lastVP.GetCBDParameters() != fVP.GetCBDParameters())
      return true;
  }

  if (lastVP.IsExplode() && (lastVP.GetExplodeFactor() != fVP.GetExplodeFactor())) return true;

  if (lastVP.IsSpecialMeshRendering()
      && (lastVP.GetSpecialMeshVolumes() != fVP.GetSpecialMeshVolumes()))
    return true;

  if (lastVP.GetTransparencyByDepth() > 0.
      && lastVP.GetTransparencyByDepthOption() != fVP.GetTransparencyByDepthOption())
    return true;

  // Time window parameters operate on the existing database so no need
  // to rebuild even if they change.

  return false;
}

void G4VulkanStoredViewer::DrawDisplayLists()
{
  // We moved these from G4VulkanViewer to G4ViewParamaters. To avoid
  // editing many lines below we introduce these convenient aliases.
#define CONVENIENT_DOUBLE_ALIAS(q) const G4double& f##q = fVP.GetTimeParameters().f##q;
#define CONVENIENT_BOOL_ALIAS(q) const G4bool& f##q = fVP.GetTimeParameters().f##q;
  CONVENIENT_DOUBLE_ALIAS(StartTime)
  CONVENIENT_DOUBLE_ALIAS(EndTime)
  CONVENIENT_DOUBLE_ALIAS(FadeFactor)
  CONVENIENT_BOOL_ALIAS(DisplayHeadTime)
  CONVENIENT_DOUBLE_ALIAS(DisplayHeadTimeX)
  CONVENIENT_DOUBLE_ALIAS(DisplayHeadTimeY)
  CONVENIENT_DOUBLE_ALIAS(DisplayHeadTimeSize)
  CONVENIENT_DOUBLE_ALIAS(DisplayHeadTimeRed)
  CONVENIENT_DOUBLE_ALIAS(DisplayHeadTimeGreen)
  CONVENIENT_DOUBLE_ALIAS(DisplayHeadTimeBlue)
  CONVENIENT_BOOL_ALIAS(DisplayLightFront)
  CONVENIENT_DOUBLE_ALIAS(DisplayLightFrontX)
  CONVENIENT_DOUBLE_ALIAS(DisplayLightFrontY)
  CONVENIENT_DOUBLE_ALIAS(DisplayLightFrontZ)
  CONVENIENT_DOUBLE_ALIAS(DisplayLightFrontT)
  CONVENIENT_DOUBLE_ALIAS(DisplayLightFrontRed)
  CONVENIENT_DOUBLE_ALIAS(DisplayLightFrontGreen)
  CONVENIENT_DOUBLE_ALIAS(DisplayLightFrontBlue)

  const G4Planes& cutaways = fVP.GetCutawayPlanes();
  G4bool cutawayUnion = fVP.IsCutaway() && fVP.GetCutawayMode() == G4ViewParameters::cutawayUnion;
  const size_t nCutaways = cutawayUnion ? cutaways.size() : 1;
  G4int iPass = 1;
  G4bool secondPassForTransparencyRequested = false;
  G4bool thirdPassForNonHiddenMarkersRequested = false;
  fDepthTestEnable = true;
  vkEnable(VK_DEPTH_TEST);
  vkDepthFunc(VK_LEQUAL);
  do
  {
    for (size_t iCutaway = 0; iCutaway < nCutaways; ++iCutaway)
    {
      if (cutawayUnion)
      {
        double a[4];
        a[0] = cutaways[iCutaway].a();
        a[1] = cutaways[iCutaway].b();
        a[2] = cutaways[iCutaway].c();
        a[3] = cutaways[iCutaway].d();
        vkClipPlane(VK_CLIP_PLANE2, a);
        vkEnable(VK_CLIP_PLANE2);
      }

      G4bool isPicking = fVP.IsPicking();

      for (size_t iPO = 0; iPO < fG4VulkanStoredSceneHandler.fPOList.size(); ++iPO)
      {
        if (POSelected(iPO))
        {
          G4VulkanStoredSceneHandler::PO& po = fG4VulkanStoredSceneHandler.fPOList[iPO];
          G4Colour c = po.fColour;
          DisplayTimePOColourModification(c, iPO);
          const G4bool isTransparent = c.GetAlpha() < 1.;
          if (iPass == 1)
          {
            if (isTransparent && transparency_enabled)
            {
              secondPassForTransparencyRequested = true;
              continue;
            }
            if (po.fMarkerOrPolyline && fVP.IsMarkerNotHidden())
            {
              thirdPassForNonHiddenMarkersRequested = true;
              continue;
            }
          }
          else if (iPass == 2)
          {  // Second pass for transparency.
            if (!isTransparent)
            {
              continue;
            }
          }
          else
          {  // Third pass for non-hidden markers
            if (!po.fMarkerOrPolyline)
            {
              continue;
            }
          }
          if (isPicking) vkLoadName(po.fPickName);
          if (transparency_enabled)
          {
            vkColor4d(c.GetRed(), c.GetGreen(), c.GetBlue(), c.GetAlpha());
          }
          else
          {
            vkColor3d(c.GetRed(), c.GetGreen(), c.GetBlue());
          }
          if (po.fMarkerOrPolyline && fVP.IsMarkerNotHidden())
          {
            if (fDepthTestEnable != false)
            {
              vkDisable(VK_DEPTH_TEST);
              fDepthTestEnable = false;
            }
          }
          else
          {
            if (fDepthTestEnable != true)
            {
              vkEnable(VK_DEPTH_TEST);
              vkDepthFunc(VK_LEQUAL);
              fDepthTestEnable = true;
            }
          }
          if (po.fpG4TextPlus)
          {
            if (po.fpG4TextPlus->fProcessing2D)
            {
              vkMatrixMode(VK_PROJECTION);
              vkPushMatrix();
              vkLoadIdentity();
              g4VkOrtho(-1., 1., -1., 1., -G4VULKAN_FLT_BIG, G4VULKAN_FLT_BIG);
              vkMatrixMode(VK_MODELVIEW);
              vkPushMatrix();
              vkLoadIdentity();
              G4VulkanTransform3D vkt(po.fTransform);
              vkMultMatrixd(vkt.GetVKMatrix());
              // This text is from a PODL. We don't want to create a new PODL.
              AddPrimitiveForASingleFrame(po.fpG4TextPlus->fG4Text);
            }
            else
            {
              vkPushMatrix();
              G4VulkanTransform3D vkt(po.fTransform);
              vkMultMatrixd(vkt.GetVKMatrix());
              // This text is from a PODL. We don't want to create a new PODL.
              AddPrimitiveForASingleFrame(po.fpG4TextPlus->fG4Text);
              vkPopMatrix();
            }

            if (po.fpG4TextPlus->fProcessing2D)
            {
              vkMatrixMode(VK_PROJECTION);
              vkPopMatrix();
              vkMatrixMode(VK_MODELVIEW);
              vkPopMatrix();
            }
          }
          else
          {
            vkPushMatrix();
            G4VulkanTransform3D vkt(po.fTransform);
            vkMultMatrixd(vkt.GetVKMatrix());
            vkCallList(po.fDisplayListId);
            vkPopMatrix();
          }
        }
      }

      G4Transform3D lastMatrixTransform;
      G4bool first = true;

      for (size_t iTO = 0; iTO < fG4VulkanStoredSceneHandler.fTOList.size(); ++iTO)
      {
        if (TOSelected(iTO))
        {
          G4VulkanStoredSceneHandler::TO& to = fG4VulkanStoredSceneHandler.fTOList[iTO];
          const G4Colour& c = to.fColour;
          const G4bool isTransparent = c.GetAlpha() < 1.;
          if (iPass == 1)
          {
            if (isTransparent && transparency_enabled)
            {
              secondPassForTransparencyRequested = true;
              continue;
            }
            if (to.fMarkerOrPolyline && fVP.IsMarkerNotHidden())
            {
              thirdPassForNonHiddenMarkersRequested = true;
              continue;
            }
          }
          else if (iPass == 2)
          {  // Second pass for transparency.
            if (!isTransparent)
            {
              continue;
            }
          }
          else
          {  // Third pass for non-hidden markers
            if (!to.fMarkerOrPolyline)
            {
              continue;
            }
          }
          if (to.fMarkerOrPolyline && fVP.IsMarkerNotHidden())
          {
            if (fDepthTestEnable != false)
            {
              vkDisable(VK_DEPTH_TEST);
              fDepthTestEnable = false;
            }
          }
          else
          {
            if (fDepthTestEnable != true)
            {
              vkEnable(VK_DEPTH_TEST);
              vkDepthFunc(VK_LEQUAL);
              fDepthTestEnable = true;
            }
          }
          if (to.fEndTime >= fStartTime && to.fStartTime <= fEndTime)
          {
            if (fVP.IsPicking()) vkLoadName(to.fPickName);
            if (to.fpG4TextPlus)
            {
              if (to.fpG4TextPlus->fProcessing2D)
              {
                vkMatrixMode(VK_PROJECTION);
                vkPushMatrix();
                vkLoadIdentity();
                g4VkOrtho(-1., 1., -1., 1., -G4VULKAN_FLT_BIG, G4VULKAN_FLT_BIG);
                vkMatrixMode(VK_MODELVIEW);
                vkPushMatrix();
                vkLoadIdentity();
              }
              G4VulkanTransform3D vkt(to.fTransform);
              vkMultMatrixd(vkt.GetVKMatrix());
              // This text is from a TODL. We don't want to create a new TODL.
              AddPrimitiveForASingleFrame(to.fpG4TextPlus->fG4Text);
              if (to.fpG4TextPlus->fProcessing2D)
              {
                vkMatrixMode(VK_PROJECTION);
                vkPopMatrix();
                vkMatrixMode(VK_MODELVIEW);
                vkPopMatrix();
              }
            }
            else
            {
              if (to.fTransform != lastMatrixTransform)
              {
                if (!first)
                {
                  vkPopMatrix();
                }
                first = false;
                vkPushMatrix();
                G4VulkanTransform3D vkt(to.fTransform);
                vkMultMatrixd(vkt.GetVKMatrix());
              }
              const G4Colour& cc = to.fColour;
              if (fFadeFactor > 0. && to.fEndTime < fEndTime)
              {
                // Brightness scaling factor
                G4double bsf =
                  1. - fFadeFactor * ((fEndTime - to.fEndTime) / (fEndTime - fStartTime));
                const G4Colour& bg = fVP.GetBackgroundColour();
                if (transparency_enabled)
                {
                  vkColor4d(bsf * cc.GetRed() + (1. - bsf) * bg.GetRed(),
                            bsf * cc.GetGreen() + (1. - bsf) * bg.GetGreen(),
                            bsf * cc.GetBlue() + (1. - bsf) * bg.GetBlue(),
                            bsf * cc.GetAlpha() + (1. - bsf) * bg.GetAlpha());
                }
                else
                {
                  vkColor3d(bsf * cc.GetRed() + (1. - bsf) * bg.GetRed(),
                            bsf * cc.GetGreen() + (1. - bsf) * bg.GetGreen(),
                            bsf * cc.GetBlue() + (1. - bsf) * bg.GetBlue());
                }
              }
              else
              {
                if (transparency_enabled)
                {
                  vkColor4d(cc.GetRed(), cc.GetGreen(), cc.GetBlue(), cc.GetAlpha());
                }
                else
                {
                  vkColor3d(cc.GetRed(), cc.GetGreen(), cc.GetBlue());
                }
              }
              vkCallList(to.fDisplayListId);
            }
            if (to.fTransform != lastMatrixTransform)
            {
              lastMatrixTransform = to.fTransform;
            }
          }
        }
      }
      if (!first)
      {
        vkPopMatrix();
      }

      if (cutawayUnion) vkDisable(VK_CLIP_PLANE2);
    }  // iCutaway

    if (iPass == 2) secondPassForTransparencyRequested = false;  // Done.
    if (iPass == 3) thirdPassForNonHiddenMarkersRequested = false;  // Done.

    if (secondPassForTransparencyRequested)
      iPass = 2;
    else if (thirdPassForNonHiddenMarkersRequested)
      iPass = 3;
    else
      break;

  } while (true);

  // Display time at "head" of time range, which is fEndTime...
  if (fDisplayHeadTime && fEndTime < G4VisAttributes::fVeryLongTime)
  {
    vkMatrixMode(VK_PROJECTION);
    vkPushMatrix();
    vkLoadIdentity();
    g4VkOrtho(-1., 1., -1., 1., -G4VULKAN_FLT_BIG, G4VULKAN_FLT_BIG);
    vkMatrixMode(VK_MODELVIEW);
    vkPushMatrix();
    vkLoadIdentity();
    G4Text headTimeText(G4BestUnit(fEndTime, "Time"),
                        G4Point3D(fDisplayHeadTimeX, fDisplayHeadTimeY, 0.));
    headTimeText.SetScreenSize(fDisplayHeadTimeSize);
    G4VisAttributes visAtts(
      G4Colour(fDisplayHeadTimeRed, fDisplayHeadTimeGreen, fDisplayHeadTimeBlue));
    headTimeText.SetVisAttributes(&visAtts);
    AddPrimitiveForASingleFrame(headTimeText);
    vkMatrixMode(VK_PROJECTION);
    vkPopMatrix();
    vkMatrixMode(VK_MODELVIEW);
    vkPopMatrix();
  }

  // Display light front...
  if (fDisplayLightFront && fEndTime < G4VisAttributes::fVeryLongTime)
  {
    G4double lightFrontRadius = (fEndTime - fDisplayLightFrontT) * c_light;
    if (lightFrontRadius > 0.)
    {
      G4Point3D lightFrontCentre(fDisplayLightFrontX, fDisplayLightFrontY, fDisplayLightFrontZ);
      G4Point3D circleCentre = lightFrontCentre;
      G4double circleRadius = lightFrontRadius;
      if (fVP.GetFieldHalfAngle() > 0.)
      {
        // Perspective view.  Find horizon centre and radius...
        G4Point3D targetPoint =
          fSceneHandler.GetScene()->GetStandardTargetPoint() + fVP.GetCurrentTargetPoint();
        G4double sceneRadius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
        if (sceneRadius <= 0.) sceneRadius = 1.;
        G4double cameraDistance = fVP.GetCameraDistance(sceneRadius);
        G4Point3D cameraPosition =
          targetPoint + cameraDistance * fVP.GetViewpointDirection().unit();
        G4Vector3D lightFrontToCameraDirection = cameraPosition - lightFrontCentre;
        G4double lightFrontCentreDistance = lightFrontToCameraDirection.mag();
        /*
         G4cout << "cameraPosition: " << cameraPosition
         << ", lightFrontCentre: " << lightFrontCentre
         << ", lightFrontRadius: " << lightFrontRadius
         << ", lightFrontCentreDistance: " << lightFrontCentreDistance
         << ", dot: " << lightFrontToCameraDirection * fVP.GetViewpointDirection()
         << G4endl;
         */
        if (lightFrontToCameraDirection * fVP.GetViewpointDirection() > 0.
            && lightFrontRadius < lightFrontCentreDistance)
        {
          // Light front in front of camera...
          G4double sineHorizonAngle = lightFrontRadius / lightFrontCentreDistance;
          circleCentre =
            lightFrontCentre
            + (lightFrontRadius * sineHorizonAngle) * lightFrontToCameraDirection.unit();
          circleRadius = lightFrontRadius * std::sqrt(1. - std::pow(sineHorizonAngle, 2));
          /*
           G4cout << "sineHorizonAngle: " << sineHorizonAngle
           << ", circleCentre: " << circleCentre
           << ", circleRadius: " << circleRadius
           << G4endl;
           */
        }
        else
        {
          circleRadius = -1.;
        }
      }
      if (circleRadius > 0.)
      {
        G4Circle lightFront(circleCentre);
        lightFront.SetWorldRadius(circleRadius);
        G4VisAttributes visAtts(
          G4Colour(fDisplayLightFrontRed, fDisplayLightFrontGreen, fDisplayLightFrontBlue));
        lightFront.SetVisAttributes(visAtts);
        AddPrimitiveForASingleFrame(lightFront);
      }
    }
  }
}

void G4VulkanStoredViewer::AddPrimitiveForASingleFrame(const G4Text& text)
{
  // We don't want this to get into a display list or a TODL or a PODL so
  // use the fMemoryForDisplayLists flag.
  fG4VulkanStoredSceneHandler.fDoNotUseDisplayList = true;
  fG4VulkanStoredSceneHandler.G4VulkanStoredSceneHandler::AddPrimitive(text);
  fG4VulkanStoredSceneHandler.fDoNotUseDisplayList = false;
}

void G4VulkanStoredViewer::AddPrimitiveForASingleFrame(const G4Circle& circle)
{
  // We don't want this to get into a display list or a TODL or a PODL so
  // use the fMemoryForDisplayLists flag.
  fG4VulkanStoredSceneHandler.fDoNotUseDisplayList = true;
  fG4VulkanStoredSceneHandler.G4VulkanStoredSceneHandler::AddPrimitive(circle);
  fG4VulkanStoredSceneHandler.fDoNotUseDisplayList = false;
}
