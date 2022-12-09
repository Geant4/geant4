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
// Andrew Walkden  7th February 1997
// Class G4OpenGLStoredViewer : Encapsulates the `storedness' of
//                            an OpenGL view, for inheritance by
//                            derived (X, Xm...) classes.

#include "G4OpenGLStoredViewer.hh"

#include "G4PhysicalConstants.hh"
#include "G4OpenGLStoredSceneHandler.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4UnitsTable.hh"
#include "G4Scene.hh"
#include "G4OpenGLTransform3D.hh"

G4OpenGLStoredViewer::G4OpenGLStoredViewer
(G4OpenGLStoredSceneHandler& sceneHandler):
G4VViewer (sceneHandler, -1),
G4OpenGLViewer (sceneHandler),
fG4OpenGLStoredSceneHandler (sceneHandler),
fDepthTestEnable(true)
{
  fLastVP = fDefaultVP;  // Update in sub-class after KernelVisitDecision
}

G4OpenGLStoredViewer::~G4OpenGLStoredViewer () {}

void G4OpenGLStoredViewer::KernelVisitDecision () {
  
  // If there's a significant difference with the last view parameters
  // of either the scene handler or this viewer, trigger a rebuild.
  
  if (!fG4OpenGLStoredSceneHandler.fTopPODL ||
      CompareForKernelVisit(fLastVP)) {
    NeedKernelVisit ();
  }
}

G4bool G4OpenGLStoredViewer::CompareForKernelVisit(G4ViewParameters& lastVP) {
  
  if (
      (lastVP.GetDrawingStyle ()    != fVP.GetDrawingStyle ())    ||
      (lastVP.GetNumberOfCloudPoints()  != fVP.GetNumberOfCloudPoints())  ||
      (lastVP.IsAuxEdgeVisible ()   != fVP.IsAuxEdgeVisible ())   ||
      (lastVP.IsCulling ()          != fVP.IsCulling ())          ||
      (lastVP.IsCullingInvisible () != fVP.IsCullingInvisible ()) ||
      (lastVP.IsDensityCulling ()   != fVP.IsDensityCulling ())   ||
      (lastVP.IsCullingCovered ()   != fVP.IsCullingCovered ())   ||
      (lastVP.GetCBDAlgorithmNumber() !=
       fVP.GetCBDAlgorithmNumber())                               ||
      // Note: Section and Cutaway can reveal back-facing faces. If
      // backface culling is implemented, the image can look strange because
      // the back-facing faces are not there. For the moment, we have disabled
      // (commented out) backface culling (it seems not to affect performance -
      // in fact, performance seems to improve), so there is no problem.
      (lastVP.IsSection ()          != fVP.IsSection ())          ||
      // Section (DCUT) is NOT implemented locally so we need to visit the kernel.
      // (lastVP.IsCutaway ()          != fVP.IsCutaway ())          ||
      // Cutaways are implemented locally so we do not need to visit the kernel.
      (lastVP.IsExplode ()          != fVP.IsExplode ())          ||
      (lastVP.GetNoOfSides ()       != fVP.GetNoOfSides ())       ||
      (lastVP.GetGlobalMarkerScale()    != fVP.GetGlobalMarkerScale())    ||
      (lastVP.GetGlobalLineWidthScale() != fVP.GetGlobalLineWidthScale()) ||
      (lastVP.IsMarkerNotHidden ()  != fVP.IsMarkerNotHidden ())  ||
      (lastVP.GetDefaultVisAttributes()->GetColour() !=
       fVP.GetDefaultVisAttributes()->GetColour())                ||
      (lastVP.GetDefaultTextVisAttributes()->GetColour() !=
       fVP.GetDefaultTextVisAttributes()->GetColour())            ||
      (lastVP.GetBackgroundColour ()!= fVP.GetBackgroundColour ())||
      (lastVP.IsPicking ()          != fVP.IsPicking ())          ||
      (lastVP.GetVisAttributesModifiers() !=
       fVP.GetVisAttributesModifiers())                           ||
      (lastVP.IsSpecialMeshRendering() !=
       fVP.IsSpecialMeshRendering())                              ||
      (lastVP.GetSpecialMeshRenderingOption() !=
       fVP.GetSpecialMeshRenderingOption())
      )
  return true;
  
  if (lastVP.IsDensityCulling () &&
      (lastVP.GetVisibleDensity () != fVP.GetVisibleDensity ()))
  return true;
  
//  /**************************************************************
//   If section (DCUT) is implemented locally, comment this out.
  if (lastVP.IsSection () &&
      (lastVP.GetSectionPlane () != fVP.GetSectionPlane ()))
    return true;
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
  
  if (lastVP.GetCBDAlgorithmNumber() > 0) {
    if (lastVP.GetCBDParameters().size() != fVP.GetCBDParameters().size()) return true;
    else if (lastVP.GetCBDParameters() != fVP.GetCBDParameters()) return true;
  }

  if (lastVP.IsExplode () &&
      (lastVP.GetExplodeFactor () != fVP.GetExplodeFactor ()))
    return true;

  if (lastVP.IsSpecialMeshRendering() &&
      (lastVP.GetSpecialMeshVolumes() != fVP.GetSpecialMeshVolumes()))
    return true;

  // Time window parameters operate on the existing database so no need
  // to rebuild even if they change.
  
  return false;
}

void G4OpenGLStoredViewer::DrawDisplayLists () {
  
  // We moved these from G4OpenGLViewer to G4ViewParamaters. To avoid
  // editing many lines below we introduce these convenient aliases.
#define CONVENIENT_DOUBLE_ALIAS(q) const G4double& f##q = fVP.Get##q();
#define CONVENIENT_BOOL_ALIAS(q) const G4bool& f##q = fVP.Is##q();
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
  G4bool cutawayUnion = fVP.IsCutaway() &&
  fVP.GetCutawayMode() == G4ViewParameters::cutawayUnion;
  const size_t nCutaways = cutawayUnion? cutaways.size(): 1;
  G4int iPass = 1;
  G4bool secondPassForTransparencyRequested = false;
  G4bool thirdPassForNonHiddenMarkersRequested = false;
  fDepthTestEnable = true;
  glEnable (GL_DEPTH_TEST); glDepthFunc (GL_LEQUAL);
  do {
    for (size_t iCutaway = 0; iCutaway < nCutaways; ++iCutaway) {
      
      if (cutawayUnion) {
        double a[4];
        a[0] = cutaways[iCutaway].a();
        a[1] = cutaways[iCutaway].b();
        a[2] = cutaways[iCutaway].c();
        a[3] = cutaways[iCutaway].d();
        glClipPlane (GL_CLIP_PLANE2, a);
        glEnable (GL_CLIP_PLANE2);
      }
      
      G4bool isPicking = fVP.IsPicking();
      
      for (size_t iPO = 0;
           iPO < fG4OpenGLStoredSceneHandler.fPOList.size(); ++iPO) {
        if (POSelected(iPO)) {
          G4OpenGLStoredSceneHandler::PO& po =
          fG4OpenGLStoredSceneHandler.fPOList[iPO];
          G4Colour c = po.fColour;
          DisplayTimePOColourModification(c,iPO);
          const G4bool isTransparent = c.GetAlpha() < 1.;
          if ( iPass == 1) {
            if (isTransparent && transparency_enabled) {
              secondPassForTransparencyRequested = true;
              continue;
            }
            if (po.fMarkerOrPolyline && fVP.IsMarkerNotHidden()) {
              thirdPassForNonHiddenMarkersRequested = true;
              continue;
            }
          } else if (iPass == 2) {  // Second pass for transparency.
            if (!isTransparent) {
              continue;
            }
          } else {  // Third pass for non-hidden markers
            if (!po.fMarkerOrPolyline) {
              continue;
            }
          }
          if (isPicking) glLoadName(po.fPickName);
          if (transparency_enabled) {
            glColor4d(c.GetRed(),c.GetGreen(),c.GetBlue(),c.GetAlpha());
          } else {
            glColor3d(c.GetRed(),c.GetGreen(),c.GetBlue());
          }
          if (po.fMarkerOrPolyline && fVP.IsMarkerNotHidden()) {
            if (fDepthTestEnable !=false) {
              glDisable (GL_DEPTH_TEST);
              fDepthTestEnable = false;
            }
          } else {
            if (fDepthTestEnable !=true) {
              glEnable (GL_DEPTH_TEST); glDepthFunc (GL_LEQUAL);
              fDepthTestEnable = true;
            }
          }
          if (po.fpG4TextPlus) {
            if (po.fpG4TextPlus->fProcessing2D) {
              glMatrixMode (GL_PROJECTION);
              glPushMatrix();
              glLoadIdentity();
              g4GlOrtho (-1., 1., -1., 1., -G4OPENGL_FLT_BIG, G4OPENGL_FLT_BIG);
              glMatrixMode (GL_MODELVIEW);
              glPushMatrix();
              glLoadIdentity();
              G4OpenGLTransform3D oglt (po.fTransform);
              glMultMatrixd (oglt.GetGLMatrix ());
              // This text is from a PODL. We don't want to create a new PODL.
              AddPrimitiveForASingleFrame(po.fpG4TextPlus->fG4Text);
            } else {
              glPushMatrix();
              G4OpenGLTransform3D oglt (po.fTransform);
              glMultMatrixd (oglt.GetGLMatrix ());
              // This text is from a PODL. We don't want to create a new PODL.
              AddPrimitiveForASingleFrame(po.fpG4TextPlus->fG4Text);
              glPopMatrix();
            }
            
            if (po.fpG4TextPlus->fProcessing2D) {
              glMatrixMode (GL_PROJECTION);
              glPopMatrix();
              glMatrixMode (GL_MODELVIEW);
              glPopMatrix();
            }
          } else {
            glPushMatrix();
            G4OpenGLTransform3D oglt (po.fTransform);
            glMultMatrixd (oglt.GetGLMatrix ());
            glCallList (po.fDisplayListId);
            glPopMatrix();
          }
        }
      }
      
      G4Transform3D lastMatrixTransform;
      G4bool first = true;

      for (size_t iTO = 0;
           iTO < fG4OpenGLStoredSceneHandler.fTOList.size(); ++iTO) {
        if (TOSelected(iTO)) {
          G4OpenGLStoredSceneHandler::TO& to =
          fG4OpenGLStoredSceneHandler.fTOList[iTO];
          const G4Colour& c = to.fColour;
          const G4bool isTransparent = c.GetAlpha() < 1.;
          if ( iPass == 1) {
            if (isTransparent && transparency_enabled) {
              secondPassForTransparencyRequested = true;
              continue;
            }
            if (to.fMarkerOrPolyline && fVP.IsMarkerNotHidden()) {
              thirdPassForNonHiddenMarkersRequested = true;
              continue;
            }
          } else if (iPass == 2) {  // Second pass for transparency.
            if (!isTransparent) {
              continue;
            }
          } else {  // Third pass for non-hidden markers
            if (!to.fMarkerOrPolyline) {
              continue;
            }
          }
          if (to.fMarkerOrPolyline && fVP.IsMarkerNotHidden()) {
            if (fDepthTestEnable !=false) {
              glDisable (GL_DEPTH_TEST);
              fDepthTestEnable = false;
            }
          } else {
            if (fDepthTestEnable !=true) {
              glEnable (GL_DEPTH_TEST); glDepthFunc (GL_LEQUAL);
              fDepthTestEnable = true;
            }
          }
          if (to.fEndTime >= fStartTime && to.fStartTime <= fEndTime) {
            if (fVP.IsPicking()) glLoadName(to.fPickName);
            if (to.fpG4TextPlus) {
              if (to.fpG4TextPlus->fProcessing2D) {
                glMatrixMode (GL_PROJECTION);
                glPushMatrix();
                glLoadIdentity();
                g4GlOrtho (-1., 1., -1., 1., -G4OPENGL_FLT_BIG, G4OPENGL_FLT_BIG);
                glMatrixMode (GL_MODELVIEW);
                glPushMatrix();
                glLoadIdentity();
              }
              G4OpenGLTransform3D oglt (to.fTransform);
              glMultMatrixd (oglt.GetGLMatrix ());
              // This text is from a TODL. We don't want to create a new TODL.
              AddPrimitiveForASingleFrame(to.fpG4TextPlus->fG4Text);
              if (to.fpG4TextPlus->fProcessing2D) {
                glMatrixMode (GL_PROJECTION);
                glPopMatrix();
                glMatrixMode (GL_MODELVIEW);
                glPopMatrix();
              }
            } else {
              if (to.fTransform != lastMatrixTransform) {
                if (! first) {
                  glPopMatrix();
                }
                first = false;
                glPushMatrix();
                G4OpenGLTransform3D oglt (to.fTransform);
                glMultMatrixd (oglt.GetGLMatrix ());
              }
              const G4Colour& cc = to.fColour;
              if (fFadeFactor > 0. && to.fEndTime < fEndTime) {
                // Brightness scaling factor
                G4double bsf = 1. - fFadeFactor *
                ((fEndTime - to.fEndTime) / (fEndTime - fStartTime));
                const G4Colour& bg = fVP.GetBackgroundColour();
                if (transparency_enabled) {
                  glColor4d
                  (bsf * cc.GetRed() + (1. - bsf) * bg.GetRed(),
                   bsf * cc.GetGreen() + (1. - bsf) * bg.GetGreen(),
                   bsf * cc.GetBlue() + (1. - bsf) * bg.GetBlue(),
                   bsf * cc.GetAlpha() + (1. - bsf) * bg.GetAlpha());
                } else {
                  glColor3d
                  (bsf * cc.GetRed() + (1. - bsf) * bg.GetRed(),
                   bsf * cc.GetGreen() + (1. - bsf) * bg.GetGreen(),
                   bsf * cc.GetBlue() + (1. - bsf) * bg.GetBlue());
                }
              } else {
                if (transparency_enabled) {
                  glColor4d(cc.GetRed(),cc.GetGreen(),cc.GetBlue(),cc.GetAlpha());
                } else {
                  glColor3d(cc.GetRed(),cc.GetGreen(),cc.GetBlue());
                }
              }
              glCallList (to.fDisplayListId);
            }
            if (to.fTransform != lastMatrixTransform) {
              lastMatrixTransform = to.fTransform;
            }
          }
        }
      }
      if (! first) {
        glPopMatrix();
      }
      
      if (cutawayUnion) glDisable (GL_CLIP_PLANE2);
    }  // iCutaway
    
    if (iPass == 2) secondPassForTransparencyRequested = false;  // Done.
    if (iPass == 3) thirdPassForNonHiddenMarkersRequested = false;  // Done.
    
    if (secondPassForTransparencyRequested) iPass = 2;
    else if (thirdPassForNonHiddenMarkersRequested) iPass = 3;
    else break;
    
  } while (true);
  
  // Display time at "head" of time range, which is fEndTime...
  if (fDisplayHeadTime && fEndTime < G4VisAttributes::fVeryLongTime) {
    glMatrixMode (GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    g4GlOrtho (-1., 1., -1., 1., -G4OPENGL_FLT_BIG, G4OPENGL_FLT_BIG);
    glMatrixMode (GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    G4Text headTimeText(G4BestUnit(fEndTime,"Time"),
                        G4Point3D(fDisplayHeadTimeX, fDisplayHeadTimeY, 0.));
    headTimeText.SetScreenSize(fDisplayHeadTimeSize);
    G4VisAttributes visAtts (G4Colour
                             (fDisplayHeadTimeRed,
                              fDisplayHeadTimeGreen,
                              fDisplayHeadTimeBlue));
    headTimeText.SetVisAttributes(&visAtts);
    AddPrimitiveForASingleFrame(headTimeText);
    glMatrixMode (GL_PROJECTION);
    glPopMatrix();
    glMatrixMode (GL_MODELVIEW);
    glPopMatrix();
  }
  
  // Display light front...
  if (fDisplayLightFront && fEndTime < G4VisAttributes::fVeryLongTime) {
    G4double lightFrontRadius = (fEndTime - fDisplayLightFrontT) * c_light;
    if (lightFrontRadius > 0.) {
      G4Point3D lightFrontCentre(fDisplayLightFrontX, fDisplayLightFrontY, fDisplayLightFrontZ);
      G4Point3D circleCentre = lightFrontCentre;
      G4double circleRadius = lightFrontRadius;
      if (fVP.GetFieldHalfAngle() > 0.) {
        // Perspective view.  Find horizon centre and radius...
        G4Point3D targetPoint = fSceneHandler.GetScene()->GetStandardTargetPoint() +
        fVP.GetCurrentTargetPoint();
        G4double sceneRadius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
        if(sceneRadius <= 0.) sceneRadius = 1.;
        G4double cameraDistance = fVP.GetCameraDistance(sceneRadius);
        G4Point3D cameraPosition = targetPoint + cameraDistance * fVP.GetViewpointDirection().unit();
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
        if (lightFrontToCameraDirection * fVP.GetViewpointDirection() > 0. && lightFrontRadius < lightFrontCentreDistance) {
          // Light front in front of camera...
          G4double sineHorizonAngle = lightFrontRadius / lightFrontCentreDistance;
          circleCentre = lightFrontCentre + (lightFrontRadius * sineHorizonAngle) * lightFrontToCameraDirection.unit();
          circleRadius = lightFrontRadius * std::sqrt(1. - std::pow(sineHorizonAngle, 2));
          /*
           G4cout << "sineHorizonAngle: " << sineHorizonAngle
           << ", circleCentre: " << circleCentre
           << ", circleRadius: " << circleRadius
           << G4endl;
           */
        } else {
          circleRadius = -1.;
        }
      }
      if (circleRadius > 0.) {
        G4Circle lightFront(circleCentre);
        lightFront.SetWorldRadius(circleRadius);
        G4VisAttributes visAtts(G4Colour
         (fDisplayLightFrontRed,
          fDisplayLightFrontGreen,
          fDisplayLightFrontBlue));
        lightFront.SetVisAttributes(visAtts);
        AddPrimitiveForASingleFrame(lightFront);
      }
    }
  }
}

void G4OpenGLStoredViewer::AddPrimitiveForASingleFrame(const G4Text& text)
{
  // We don't want this to get into a display list or a TODL or a PODL so
  // use the fMemoryForDisplayLists flag.
  fG4OpenGLStoredSceneHandler.fDoNotUseDisplayList = true;
  fG4OpenGLStoredSceneHandler.G4OpenGLStoredSceneHandler::AddPrimitive(text);
  fG4OpenGLStoredSceneHandler.fDoNotUseDisplayList = false;
}

void G4OpenGLStoredViewer::AddPrimitiveForASingleFrame(const G4Circle& circle)
{
  // We don't want this to get into a display list or a TODL or a PODL so
  // use the fMemoryForDisplayLists flag.
  fG4OpenGLStoredSceneHandler.fDoNotUseDisplayList = true;
  fG4OpenGLStoredSceneHandler.G4OpenGLStoredSceneHandler::AddPrimitive(circle);
  fG4OpenGLStoredSceneHandler.fDoNotUseDisplayList = false;
}
