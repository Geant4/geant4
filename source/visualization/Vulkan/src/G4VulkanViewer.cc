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
// Ported to Vulkan from the OpenGL work by Andrew Walkden (27th March 1996).
// Vulkan view - opens window, hard copy, etc.

#include "G4VulkanViewer.hh"

#include "G4VulkanSceneHandler.hh"
#include "G4VulkanTransform3D.hh"
#include "G4vk2ps.hh"
#define VK2PS_TEXT_B TOOLS_GL2PS_TEXT_B
#define VK2PS_TEXT_BL TOOLS_GL2PS_TEXT_BL
#define VK2PS_TEXT_BR TOOLS_GL2PS_TEXT_BR

#include "G4AttCheck.hh"
#include "G4AttHolder.hh"
#include "G4LogicalVolume.hh"
#include "G4Normal3D.hh"
#include "G4Plane3D.hh"
#include "G4Point3D.hh"
#include "G4Scene.hh"
#include "G4Text.hh"
#include "G4VSolid.hh"
#include "G4VisExtent.hh"

#include <iomanip>
#include <sstream>
#include <string>

G4VulkanViewer::G4VulkanViewer(G4VulkanSceneHandler& scene)
  : G4VViewer(scene, -1),
    fPrintColour(true),
    fVectoredPs(true),
    fVulkanSceneHandler(scene),
    background(G4Colour(0., 0., 0.)),
    transparency_enabled(true),
    antialiasing_enabled(false),
    haloing_enabled(false),
    fRot_sens(1.),
    fPan_sens(0.01),
    fWinSize_x(0),
    fWinSize_y(0),
    fDefaultExportImageFormat("pdf"),
    fExportImageFormat("pdf"),
    fExportFilenameIndex(0),
    fPrintSizeX(-1),
    fPrintSizeY(-1),
    fPointSize(0),
    fDefaultExportFilename("G4Vulkan"),
    fSizeHasChanged(0),
    fVk2psDefaultLineWith(1),
    fVk2psDefaultPointSize(2),
    fVkViewInitialized(false),
    fIsGettingPickInfos(false)
{
  // Make changes to view parameters for Vulkan...
  fVP.SetAutoRefresh(true);
  fDefaultVP.SetAutoRefresh(true);
  fVK2PSAction = new G4vk2ps();
  tools_gl2ps_gl_funcs_t _funcs = {(tools_glIsEnabled_func)vkIsEnabled,
                                   (tools_glBegin_func)vkBegin,
                                   (tools_glEnd_func)vkEnd,
                                   (tools_glGetFloatv_func)vkGetFloatv,
                                   (tools_glVertex3f_func)vkVertex3f,
                                   (tools_glGetBooleanv_func)vkGetBooleanv,
                                   (tools_glGetIntegerv_func)vkGetIntegerv,
                                   (tools_glRenderMode_func)vkRenderMode,
                                   (tools_glFeedbackBuffer_func)vkFeedbackBuffer,
                                   (tools_glPassThrough_func)vkPassThrough};
  fVK2PSAction->setVulkanFunctions(&_funcs);

  // add supported export image format
  addExportImageFormat("eps");
  addExportImageFormat("ps");
  addExportImageFormat("pdf");
  addExportImageFormat("svg");

  // Change the default name
  fExportFilename += fDefaultExportFilename + "_" + GetShortName().data();

  //  vkClearColor (0.0, 0.0, 0.0, 0.0);
  //  vkClearDepth (1.0);
  //  vkDisable (VK_BLEND);
  //  vkDisable (VK_LINE_SMOOTH);
  //  vkDisable (VK_POLYGON_SMOOTH);
}

G4VulkanViewer::~G4VulkanViewer()
{
  delete fVK2PSAction;
}

void G4VulkanViewer::InitializeVKView()
{
  if (fWinSize_x == 0)
  {
    fWinSize_x = fVP.GetWindowSizeHintX();
  }
  if (fWinSize_y == 0)
  {
    fWinSize_y = fVP.GetWindowSizeHintY();
  }

  vkClearColor(0.0, 0.0, 0.0, 0.0);
  vkClearDepth(1.0);
  vkDisable(VK_LINE_SMOOTH);
  vkDisable(VK_POLYGON_SMOOTH);
  vkDisable(VK_POINT_SMOOTH);

  // clear the buffers and window?
  ClearView();
  FinishView();

  vkDepthFunc(VK_LEQUAL);
  vkDepthMask(VK_TRUE);

  vkEnable(VK_BLEND);
  vkBlendFunc(VK_SRC_ALPHA, VK_ONE_MINUS_SRC_ALPHA);
}

void G4VulkanViewer::ClearView()
{
  ClearViewWithoutFlush();

  if (!isFramebufferReady())
  {
    return;
  }

  vkFlush();
}

void G4VulkanViewer::ClearViewWithoutFlush()
{
  // Ready for clear ?
  // See : http://lists.apple.com/archives/mac-opengl/2012/Jul/msg00038.html
  if (!isFramebufferReady())
  {
    return;
  }

  vkClearColor(background.GetRed(), background.GetGreen(), background.GetBlue(), 1.);
  vkClearDepth(1.0);
  // Below line does not compile with Mesa includes.
  // vkClear (VK_COLOR_BUFFER_BIT | VK_DEPTH_BUFFER_BIT);
  vkClear(VK_COLOR_BUFFER_BIT);
  vkClear(VK_DEPTH_BUFFER_BIT);
  vkClear(VK_STENCIL_BUFFER_BIT);
}

void G4VulkanViewer::ResizeWindow(unsigned int aWidth, unsigned int aHeight)
{
  if ((fWinSize_x != aWidth) || (fWinSize_y != aHeight))
  {
    fWinSize_x = aWidth;
    fWinSize_y = aHeight;
    fSizeHasChanged = true;
  }
  else
  {
    fSizeHasChanged = false;
  }
}

/**
 * Set the viewport of the scene
 * MAXIMUM SIZE is :
 * VKint dims[2];
 * vkGetIntegerv(VK_MAX_VIEWPORT_DIMS, dims);
 */
void G4VulkanViewer::ResizeVKView()
{
  // Check size
  VKint dims[2];
  dims[0] = 0;
  dims[1] = 0;

  vkGetIntegerv(VK_MAX_VIEWPORT_DIMS, dims);

  if ((dims[0] != 0) && (dims[1] != 0))
  {
    if (fWinSize_x > (unsigned)dims[0])
    {
      G4cerr << "Try to resize view greater than max X viewport dimension. Desired size "
             << fWinSize_x << " is resize to " << dims[0] << G4endl;
      fWinSize_x = dims[0];
    }
    if (fWinSize_y > (unsigned)dims[1])
    {
      G4cerr << "Try to resize view greater than max Y viewport dimension. Desired size "
             << fWinSize_y << " is resize to " << dims[1] << G4endl;
      fWinSize_y = dims[1];
    }
  }

  vkViewport(0, 0, fWinSize_x, fWinSize_y);
}

void G4VulkanViewer::SetView()
{
  // if getting pick infos, should not resize the view.
  if (fIsGettingPickInfos) return;

  if (!fSceneHandler.GetScene())
  {
    return;
  }
  // Calculates view representation based on extent of object being
  // viewed and (initial) viewpoint.  (Note: it can change later due
  // to user interaction via visualization system's GUI.)

  // Lighting.
  VKfloat lightPosition[4];
  lightPosition[0] = fVP.GetActualLightpointDirection().x();
  lightPosition[1] = fVP.GetActualLightpointDirection().y();
  lightPosition[2] = fVP.GetActualLightpointDirection().z();
  lightPosition[3] = 0.;
  // Light position is "true" light direction, so must come after gluLookAt.
  VKfloat ambient[] = {0.2f, 0.2f, 0.2f, 1.f};
  VKfloat diffuse[] = {0.8f, 0.8f, 0.8f, 1.f};
  vkEnable(VK_LIGHT0);
  vkLightfv(VK_LIGHT0, VK_AMBIENT, ambient);
  vkLightfv(VK_LIGHT0, VK_DIFFUSE, diffuse);

  G4double ratioX = 1;
  G4double ratioY = 1;
  if (fWinSize_y > fWinSize_x)
  {
    ratioX = ((G4double)fWinSize_y) / ((G4double)fWinSize_x);
  }
  if (fWinSize_x > fWinSize_y)
  {
    ratioY = ((G4double)fWinSize_x) / ((G4double)fWinSize_y);
  }

  // Get radius of scene, etc.
  // Note that this procedure properly takes into account zoom, dolly and pan.
  const G4Point3D targetPoint =
    fSceneHandler.GetScene()->GetStandardTargetPoint() + fVP.GetCurrentTargetPoint();
  G4double radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
  if (radius <= 0.) radius = 1.;
  const G4double cameraDistance = fVP.GetCameraDistance(radius);
  const G4Point3D cameraPosition =
    targetPoint + cameraDistance * fVP.GetViewpointDirection().unit();
  const VKdouble pnear = fVP.GetNearDistance(cameraDistance, radius);
  const VKdouble pfar = fVP.GetFarDistance(cameraDistance, pnear, radius);
  const VKdouble right = fVP.GetFrontHalfHeight(pnear, radius) * ratioY;
  const VKdouble left = -right;
  const VKdouble top = fVP.GetFrontHalfHeight(pnear, radius) * ratioX;
  const VKdouble bottom = -top;

  // FIXME
  ResizeVKView();
  // SHOULD SetWindowsSizeHint()...

  vkMatrixMode(VK_PROJECTION);  // set up Frustum.
  vkLoadIdentity();

  const G4Vector3D scaleFactor = fVP.GetScaleFactor();
  vkScaled(scaleFactor.x(), scaleFactor.y(), scaleFactor.z());

  if (fVP.GetFieldHalfAngle() == 0.)
  {
    g4VkOrtho(left, right, bottom, top, pnear, pfar);
  }
  else
  {
    g4VkFrustum(left, right, bottom, top, pnear, pfar);
  }

  vkMatrixMode(VK_MODELVIEW);  // apply further transformations to scene.
  vkLoadIdentity();

  const G4Normal3D& upVector = fVP.GetUpVector();
  G4Point3D vktarget;
  if (cameraDistance > 1.e-6 * radius)
  {
    vktarget = targetPoint;
  }
  else
  {
    vktarget = targetPoint - radius * fVP.GetViewpointDirection().unit();
  }

  const G4Point3D& pCamera = cameraPosition;  // An alias for brevity.

  g4VkLookAt(pCamera.x(), pCamera.y(), pCamera.z(),  // Viewpoint.
              vktarget.x(), vktarget.y(), vktarget.z(),  // Target point.
              upVector.x(), upVector.y(), upVector.z());  // Up vector.
  // Light position is "true" light direction, so must come after gluLookAt.
  vkLightfv(VK_LIGHT0, VK_POSITION, lightPosition);

  // The idea is to use back-to-back clipping planes.  This can cut an object
  // down to just a few pixels, which can make it difficult to see.  So, for
  // now, comment this out and use the generic (Boolean) method, via
  // G4VSolid* G4VulkanSceneHandler::CreateSectionSolid ()
  // { return G4VSceneHandler::CreateSectionSolid(); }
  //  if (fVP.IsSection () ) {  // pair of back to back clip planes.
  //    const G4Plane3D& sp = fVP.GetSectionPlane ();
  //    double sArray[4];
  //    sArray[0] = sp.a();
  //    sArray[1] = sp.b();
  //    sArray[2] = sp.c();
  //    sArray[3] = sp.d() + radius * 1.e-05;
  //    vkClipPlane (VK_CLIP_PLANE0, sArray);
  //    vkEnable (VK_CLIP_PLANE0);
  //    sArray[0] = -sp.a();
  //    sArray[1] = -sp.b();
  //    sArray[2] = -sp.c();
  //    sArray[3] = -sp.d() + radius * 1.e-05;
  //    vkClipPlane (VK_CLIP_PLANE1, sArray);
  //    vkEnable (VK_CLIP_PLANE1);
  //  } else {
  //    vkDisable (VK_CLIP_PLANE0);
  //    vkDisable (VK_CLIP_PLANE1);
  //  }

  // What we call intersection of cutaways is easy in Vulkan.  You
  // just keep cutting.  Unions are more tricky - you have to have
  // multiple passes and this is handled in
  // G4VulkanImmediate/StoredViewer::ProcessView.
  const G4Planes& cutaways = fVP.GetCutawayPlanes();
  size_t nPlanes = cutaways.size();
  if (fVP.IsCutaway() && fVP.GetCutawayMode() == G4ViewParameters::cutawayIntersection)
  {
    double a[4];
    a[0] = cutaways[0].a();
    a[1] = cutaways[0].b();
    a[2] = cutaways[0].c();
    a[3] = cutaways[0].d();
    vkClipPlane(VK_CLIP_PLANE2, a);
    vkEnable(VK_CLIP_PLANE2);
    if (nPlanes > 1)
    {
      a[0] = cutaways[1].a();
      a[1] = cutaways[1].b();
      a[2] = cutaways[1].c();
      a[3] = cutaways[1].d();
      vkClipPlane(VK_CLIP_PLANE3, a);
      vkEnable(VK_CLIP_PLANE3);
    }
    if (nPlanes > 2)
    {
      a[0] = cutaways[2].a();
      a[1] = cutaways[2].b();
      a[2] = cutaways[2].c();
      a[3] = cutaways[2].d();
      vkClipPlane(VK_CLIP_PLANE4, a);
      vkEnable(VK_CLIP_PLANE4);
    }
  }
  else
  {
    vkDisable(VK_CLIP_PLANE2);
    vkDisable(VK_CLIP_PLANE3);
    vkDisable(VK_CLIP_PLANE4);
  }

  // Background.
  background = fVP.GetBackgroundColour();
}

void G4VulkanViewer::ResetView()
{
  G4VViewer::ResetView();
  fRot_sens = 1;
  fPan_sens = 0.01;
}

void G4VulkanViewer::HaloingFirstPass()
{
  // To perform haloing, first Draw all information to the depth buffer
  // alone, using a chunky line width, and then Draw all info again, to
  // the colour buffer, setting a thinner line width an the depth testing
  // function to less than or equal, so if two lines cross, the one
  // passing behind the other will not pass the depth test, and so not
  // get rendered either side of the infront line for a short distance.

  // First, disable writing to the colo(u)r buffer...
  vkColorMask(VK_FALSE, VK_FALSE, VK_FALSE, VK_FALSE);

  // Now enable writing to the depth buffer...
  vkDepthMask(VK_TRUE);
  vkDepthFunc(VK_LESS);
  vkClearDepth(1.0);

  // Finally, set the line width to something wide...
  ChangeLineWidth(3.0);
}

void G4VulkanViewer::HaloingSecondPass()
{
  // And finally, turn the colour buffer back on with a sesible line width...
  vkColorMask(VK_TRUE, VK_TRUE, VK_TRUE, VK_TRUE);
  vkDepthFunc(VK_LEQUAL);
  ChangeLineWidth(1.0);
}

G4String G4VulkanViewer::Pick(VKdouble x, VKdouble y)
{
  const std::vector<G4VulkanViewerPickMap*>& pickMap = GetPickDetails(x, y);
  G4String txt = "";
  if (pickMap.size() == 0)
  {
    //        txt += "No hits recorded.";;
  }
  else
  {
    for (unsigned int a = 0; a < pickMap.size(); a++)
    {
      if (pickMap[a]->getAttributes().size() > 0)
      {
        txt += pickMap[a]->print();
      }
    }
  }
  return txt;
}

const std::vector<G4VulkanViewerPickMap*>& G4VulkanViewer::GetPickDetails(VKdouble x, VKdouble y)
{
  static std::vector<G4VulkanViewerPickMap*> pickMapVector;
  for (auto pickMap : pickMapVector)
  {
    delete pickMap;
  }
  pickMapVector.clear();

  const G4int BUFSIZE = 512;
  VKuint selectBuffer[BUFSIZE];
  vkSelectBuffer(BUFSIZE, selectBuffer);
  vkRenderMode(VK_SELECT);
  vkInitNames();
  vkPushName(0);
  vkMatrixMode(VK_PROJECTION);
  G4double currentProjectionMatrix[16];
  vkGetDoublev(VK_PROJECTION_MATRIX, currentProjectionMatrix);
  vkPushMatrix();
  vkLoadIdentity();
  VKint viewport[4];
  vkGetIntegerv(VK_VIEWPORT, viewport);
  /*  G4cout
    << "viewport, x,y: "
    << viewport[0] << ',' << viewport[1] << ',' << viewport[2] << ',' << viewport[3]
    << ", " << x << ',' << y
    << G4endl;
  */
  fIsGettingPickInfos = true;
  // Define 5x5 pixel pick area
  g4VkPickMatrix(x, viewport[3] - y, 5., 5., viewport);
  vkMultMatrixd(currentProjectionMatrix);
  vkMatrixMode(VK_MODELVIEW);
  DrawView();
  VKint hits = vkRenderMode(VK_RENDER);
  fIsGettingPickInfos = false;
  if (hits < 0)
  {
    G4cout << "Too many hits.  Zoom in to reduce overlaps." << G4endl;
    goto restoreMatrices;
  }
  if (hits > 0)
  {
    VKuint* p = selectBuffer;
    for (VKint i = 0; i < hits; ++i)
    {
      VKuint nnames = *p++;
      // This bit of debug code or...
      // VKuint zmin = *p++;
      // VKuint zmax = *p++;
      // G4cout << "Hit " << i << ": " << nnames << " names"
      //       << "\nzmin: " << zmin << ", zmax: " << zmax << G4endl;
      // ...just increment the pointer
      p++;
      p++;
      for (VKuint j = 0; j < nnames; ++j)
      {
        VKuint name = *p++;
        std::map<VKuint, G4AttHolder*>::iterator iter = fVulkanSceneHandler.fPickMap.find(name);
        if (iter != fVulkanSceneHandler.fPickMap.end())
        {
          G4AttHolder* attHolder = iter->second;
          if (attHolder && attHolder->GetAttDefs().size())
          {
            for (size_t iAtt = 0; iAtt < attHolder->GetAttDefs().size(); ++iAtt)
            {
              std::ostringstream oss;
              oss << G4AttCheck(attHolder->GetAttValues()[iAtt], attHolder->GetAttDefs()[iAtt]);
              G4VulkanViewerPickMap* pickMap = new G4VulkanViewerPickMap();
              //              G4cout
              //              << "i,j, attHolder->GetAttDefs().size(): "
              //              << i << ',' << j
              //              << ", " << attHolder->GetAttDefs().size()
              //              << G4endl;
              //              G4cout << "G4VulkanViewer::GetPickDetails: " << oss.str() << G4endl;
              pickMap->addAttributes(oss.str());
              pickMap->setHitNumber(i);
              pickMap->setSubHitNumber(j);
              pickMap->setPickName(name);
              pickMapVector.push_back(pickMap);
            }
          }
        }
      }
    }
  }

restoreMatrices:
  vkMatrixMode(VK_PROJECTION);
  vkPopMatrix();
  vkMatrixMode(VK_MODELVIEW);

  return pickMapVector;
}

VKubyte* G4VulkanViewer::grabPixels(int inColor, unsigned int width, unsigned int height)
{
  VKubyte* buffer;
  VKint swapbytes, lsbfirst, rowlength;
  VKint skiprows, skippixels, alignment;
  VKenum format;
  int size;

  if (inColor)
  {
    format = VK_RGB;
    size = width * height * 3;
  }
  else
  {
    format = VK_LUMINANCE;
    size = width * height * 1;
  }

  buffer = new VKubyte[size];
  if (buffer == NULL) return NULL;

  vkGetIntegerv(VK_UNPACK_SWAP_BYTES, &swapbytes);
  vkGetIntegerv(VK_UNPACK_LSB_FIRST, &lsbfirst);
  vkGetIntegerv(VK_UNPACK_ROW_LENGTH, &rowlength);

  vkGetIntegerv(VK_UNPACK_SKIP_ROWS, &skiprows);
  vkGetIntegerv(VK_UNPACK_SKIP_PIXELS, &skippixels);
  vkGetIntegerv(VK_UNPACK_ALIGNMENT, &alignment);

  vkPixelStorei(VK_UNPACK_SWAP_BYTES, VK_FALSE);
  vkPixelStorei(VK_UNPACK_LSB_FIRST, VK_FALSE);
  vkPixelStorei(VK_UNPACK_ROW_LENGTH, 0);

  vkPixelStorei(VK_UNPACK_SKIP_ROWS, 0);
  vkPixelStorei(VK_UNPACK_SKIP_PIXELS, 0);
  vkPixelStorei(VK_UNPACK_ALIGNMENT, 1);

  vkReadBuffer(VK_FRONT);
  vkReadPixels(0, 0, (VKsizei)width, (VKsizei)height, format, VK_UNSIGNED_BYTE, (VKvoid*)buffer);

  vkPixelStorei(VK_UNPACK_SWAP_BYTES, swapbytes);
  vkPixelStorei(VK_UNPACK_LSB_FIRST, lsbfirst);
  vkPixelStorei(VK_UNPACK_ROW_LENGTH, rowlength);

  vkPixelStorei(VK_UNPACK_SKIP_ROWS, skiprows);
  vkPixelStorei(VK_UNPACK_SKIP_PIXELS, skippixels);
  vkPixelStorei(VK_UNPACK_ALIGNMENT, alignment);

  return buffer;
}

bool G4VulkanViewer::printVectoredEPS()
{
  return printVk2PS();
}

bool G4VulkanViewer::printNonVectoredEPS()
{
  int width = getRealExportWidth();
  int height = getRealExportHeight();

  FILE* fp;
  VKubyte* pixels;
  VKubyte* curpix;
  int components, pos, i;

  pixels = grabPixels(fPrintColour, width, height);

  if (pixels == NULL)
  {
    G4cerr << "Failed to get pixels from Vulkan viewport" << G4endl;
    return false;
  }
  if (fPrintColour)
  {
    components = 3;
  }
  else
  {
    components = 1;
  }
  std::string name = getRealPrintFilename();
  fp = fopen(name.c_str(), "w");
  if (fp == NULL)
  {
    G4cerr << "Can't open filename " << name.c_str() << G4endl;
    return false;
  }

  fprintf(fp, "%%!PS-Adobe-2.0 EPSF-1.2\n");
  fprintf(fp, "%%%%Title: %s\n", name.c_str());
  fprintf(fp, "%%%%Creator: Vulkan pixmap render output\n");
  fprintf(fp, "%%%%BoundingBox: 0 0 %d %d\n", width, height);
  fprintf(fp, "%%%%EndComments\n");
  fprintf(fp, "gsave\n");
  fprintf(fp, "/bwproc {\n");
  fprintf(fp, "    rgbproc\n");
  fprintf(fp, "    dup length 3 idiv string 0 3 0 \n");
  fprintf(fp, "    5 -1 roll {\n");
  fprintf(fp, "    add 2 1 roll 1 sub dup 0 eq\n");
  fprintf(fp, "    { pop 3 idiv 3 -1 roll dup 4 -1 roll dup\n");
  fprintf(fp, "       3 1 roll 5 -1 roll } put 1 add 3 0 \n");
  fprintf(fp, "    { 2 1 roll } ifelse\n");
  fprintf(fp, "    }forall\n");
  fprintf(fp, "    pop pop pop\n");
  fprintf(fp, "} def\n");
  fprintf(fp, "systemdict /colorimage known not {\n");
  fprintf(fp, "   /colorimage {\n");
  fprintf(fp, "       pop\n");
  fprintf(fp, "       pop\n");
  fprintf(fp, "       /rgbproc exch def\n");
  fprintf(fp, "       { bwproc } image\n");
  fprintf(fp, "   }  def\n");
  fprintf(fp, "} if\n");
  fprintf(fp, "/picstr %d string def\n", width * components);
  fprintf(fp, "%d %d scale\n", width, height);
  fprintf(fp, "%d %d %d\n", width, height, 8);
  fprintf(fp, "[%d 0 0 %d 0 0]\n", width, height);
  fprintf(fp, "{currentfile picstr readhexstring pop}\n");
  fprintf(fp, "false %d\n", components);
  fprintf(fp, "colorimage\n");

  curpix = (VKubyte*)pixels;
  pos = 0;
  for (i = width * height * components; i > 0; i--)
  {
    fprintf(fp, "%02hx ", (unsigned short)(*(curpix++)));
    if (++pos >= 32)
    {
      fprintf(fp, "\n");
      pos = 0;
    }
  }
  if (pos) fprintf(fp, "\n");

  fprintf(fp, "grestore\n");
  fprintf(fp, "showpage\n");
  delete[] pixels;
  fclose(fp);

  // Reset for next time (useful if size change)
  //  fPrintSizeX = -1;
  //  fPrintSizeY = -1;

  return true;
}

/** Return if gl2ps is currently writing
 */
bool G4VulkanViewer::isVk2psWriting()
{
  if (!fVK2PSAction) return false;
  if (fVK2PSAction->fileWritingEnabled())
  {
    return true;
  }
  return false;
}

G4bool G4VulkanViewer::isFramebufferReady()
{
  bool check = false;
#ifdef G4VIS_BUILD_VULKANQT_DRIVER
  check = true;
#endif
#ifdef G4VIS_BUILD_VULKANX_DRIVER
  check = false;
#endif
#ifdef G4VIS_BUILD_VULKANXM_DRIVER
  check = false;
#endif
#ifdef G4VIS_BUILD_VULKANWIN32_DRIVER
  check = false;
#endif

#if VK_ARB_framebuffer_object
  if (check)
  {
    //    if ( vkCheckFramebufferStatus(VK_FRAMEBUFFER) == VK_FRAMEBUFFER_UNDEFINED) {
    //      return false;
    //    }
  }
#endif
  return true;
}

/* Draw Gl2Ps text if needed
 */
void G4VulkanViewer::DrawText(const G4Text& g4text)
{
  // gl2ps or GL window ?
  if (isVk2psWriting())
  {
    G4VSceneHandler::MarkerSizeType sizeType;
    G4double size = fSceneHandler.GetMarkerSize(g4text, sizeType);
    G4Point3D position = g4text.GetPosition();

    G4String textString = g4text.GetText();

    vkRasterPos3d(position.x(), position.y(), position.z());
    VKint align = VK2PS_TEXT_B;

    switch (g4text.GetLayout())
    {
      case G4Text::left:
        align = VK2PS_TEXT_BL;
        break;
      case G4Text::centre:
        align = VK2PS_TEXT_B;
        break;
      case G4Text::right:
        align = VK2PS_TEXT_BR;
    }

    fVK2PSAction->addTextOpt(textString.c_str(), "Times-Roman", VKshort(size), align, 0);
  }
  else
  {
    static G4int callCount = 0;
    ++callCount;
    // if (callCount <= 10 || callCount%100 == 0) {
    if (callCount <= 1)
    {
      G4cout << "G4VulkanViewer::DrawText: Not implemented for \"" << fName << "\"\n  Called with "
             << g4text << G4endl;
    }
  }
}

/** Change PointSize on gl2ps if needed
 */
void G4VulkanViewer::ChangePointSize(G4double size)
{
  if (isVk2psWriting())
  {
    fVK2PSAction->setPointSize(int(size));
  }
  else
  {
    vkPointSize(size);
  }
}

/** Change LineSize on gl2ps if needed
 */
void G4VulkanViewer::ChangeLineWidth(G4double width)
{
  if (isVk2psWriting())
  {
    fVK2PSAction->setLineWidth(int(width));
  }
  else
  {
    vkLineWidth(width);
  }
}

/**
 Export image with the given name with width and height
 Several cases :
 If name is "", filename will have the default value
 If name is "toto.png", set the name to "toto" and the format to "png". No incremented suffix is
 added. If name is "toto", set the name to "toto" and the format to default (or current format if
 specify). Will also add an incremented suffix at the end of the file
*/
bool G4VulkanViewer::exportImage(std::string name, int width, int height)
{
  if (!setExportFilename(name))
  {
    return false;
  }

  if ((width != -1) && (height != -1))
  {
    setExportSize(width, height);
  }

  if (fExportImageFormat == "eps")
  {
    fVK2PSAction->setExportImageFormat_EPS();
  }
  else if (fExportImageFormat == "ps")
  {
    fVK2PSAction->setExportImageFormat_PS();
  }
  else if (fExportImageFormat == "svg")
  {
    fVK2PSAction->setExportImageFormat_SVG();
  }
  else if (fExportImageFormat == "pdf")
  {
    fVK2PSAction->setExportImageFormat_PDF();
  }
  else
  {
    setExportImageFormat(
      fExportImageFormat,
      true);  // will display a message if this format is not correct for the current viewer
    return false;
  }

  bool res;

  // Change the LC_NUMERIC value in order to have "." separtor and not ","
  // This case is only useful for French, Canadien...
  size_t len = strlen(setlocale(LC_NUMERIC, NULL));
  char* oldLocale = (char*)(malloc(len + 1));
  if (oldLocale != NULL) strncpy(oldLocale, setlocale(LC_NUMERIC, NULL), len);
  setlocale(LC_NUMERIC, "C");

  if (((fExportImageFormat == "eps") || (fExportImageFormat == "ps")) && (!fVectoredPs))
  {
    res = printNonVectoredEPS();
  }
  else
  {
    res = printVectoredEPS();
  }

  // restore the local
  if (oldLocale)
  {
    setlocale(LC_NUMERIC, oldLocale);
    free(oldLocale);
  }

  if (res == false)
  {
    G4cerr << "Error saving file... " << getRealPrintFilename().c_str() << G4endl;
  }
  else
  {
    G4cout << "File " << getRealPrintFilename().c_str() << " size: " << getRealExportWidth() << "x"
           << getRealExportHeight() << " has been saved " << G4endl;

    // increment index if necessary
    if (fExportFilenameIndex != -1)
    {
      fExportFilenameIndex++;
    }
  }

  return res;
}

bool G4VulkanViewer::printVk2PS()
{
  int width = getRealExportWidth();
  int height = getRealExportHeight();
  bool res = true;

  // no need to redraw at each new primitive for printvk2PS
  G4VulkanSceneHandler& vkSceneHandler = dynamic_cast<G4VulkanSceneHandler&>(fSceneHandler);
  G4VulkanSceneHandler::FlushAction originalFlushAction = vkSceneHandler.GetFlushAction();
  vkSceneHandler.SetFlushAction(G4VulkanSceneHandler::never);

  if (!fVK2PSAction) return false;

  fVK2PSAction->setFileName(getRealPrintFilename().c_str());
  // try to resize
  int X = fWinSize_x;
  int Y = fWinSize_y;

  fWinSize_x = width;
  fWinSize_y = height;
  // Laurent G. 16/03/10 : Not the good way to do.
  // We should draw in a new offscreen context instead of
  // resizing and drawing in current window...
  // This should be solve when we will do an offscreen method
  // to render Vulkan
  // See :
  // http://developer.apple.com/Mac/library/documentation/GraphicsImaging/Conceptual/Vulkan-MacProgGuide/vulkan_offscreen/vulkan_offscreen.html
  // http://www.songho.ca/opengl/vk_fbo.html

  ResizeVKView();
  bool extendBuffer = true;
  bool endWriteAction = false;
  bool beginWriteAction = true;
  bool filePointerOk = true;
  while ((extendBuffer) && (!endWriteAction) && (filePointerOk))
  {
    beginWriteAction = fVK2PSAction->enableFileWriting();
    if (beginWriteAction)
    {
      VKint vp[4];
      ::vkGetIntegerv(VK_VIEWPORT, vp);
      fVK2PSAction->setViewport(vp[0], vp[1], vp[2], vp[3]);
      beginWriteAction = fVK2PSAction->beginPage();
    }

    // 3 cases :
    // - true
    // - false && ! fVK2PSAction->fileWritingEnabled() => bad file name
    // - false && fVK2PSAction->fileWritingEnabled() => buffer size problem ?

    filePointerOk = fVK2PSAction->fileWritingEnabled();

    if (beginWriteAction)
    {
      // Set the viewport
      // By default, we choose the line width (trajectories...)
      fVK2PSAction->setLineWidth(fVk2psDefaultLineWith);
      // By default, we choose the point size (markers...)
      fVK2PSAction->setPointSize(fVk2psDefaultPointSize);

      DrawView();

      endWriteAction = fVK2PSAction->endPage();
      fVK2PSAction->disableFileWriting();
    }
    if (filePointerOk)
    {
      if ((!endWriteAction) || (!beginWriteAction))
      {
        extendBuffer = fVK2PSAction->extendBufferSize();
      }
    }
  }
  fVK2PSAction->resetBufferSizeParameters();

  if (!extendBuffer)
  {
    G4cerr << "ERROR: gl2ps buffer size is not big enough to print this geometry. Try to extend "
              "it. No output produced"
           << G4endl;
    res = false;
  }
  if (!beginWriteAction)
  {
    G4cerr << "ERROR: saving file " << getRealPrintFilename().c_str()
           << ". Check read/write access. No output produced" << G4endl;
    res = false;
  }
  if (!endWriteAction)
  {
    G4cerr << "gl2ps error. No output produced" << G4endl;
    res = false;
  }
  fWinSize_x = X;
  fWinSize_y = Y;

  vkSceneHandler.SetFlushAction(originalFlushAction);

  // Reset for next time (useful is size change)
  //  fPrintSizeX = 0;
  //  fPrintSizeY = 0;

  return res;
}

unsigned int G4VulkanViewer::getWinWidth() const
{
  return fWinSize_x;
}

unsigned int G4VulkanViewer::getWinHeight() const
{
  return fWinSize_y;
}

G4bool G4VulkanViewer::sizeHasChanged()
{
  return fSizeHasChanged;
}

G4int G4VulkanViewer::getRealExportWidth()
{
  if (fPrintSizeX == -1)
  {
    return fWinSize_x;
  }
  VKint dims[2];
  vkGetIntegerv(VK_MAX_VIEWPORT_DIMS, dims);

  // L.Garnier 01-2010: Some problems with mac 10.6
  if ((dims[0] != 0) && (dims[1] != 0))
  {
    if (fPrintSizeX > dims[0])
    {
      return dims[0];
    }
  }
  if (fPrintSizeX < -1)
  {
    return 0;
  }
  return fPrintSizeX;
}

G4int G4VulkanViewer::getRealExportHeight()
{
  if (fPrintSizeY == -1)
  {
    return fWinSize_y;
  }
  VKint dims[2];
  vkGetIntegerv(VK_MAX_VIEWPORT_DIMS, dims);

  // L.Garnier 01-2010: Some problems with mac 10.6
  if ((dims[0] != 0) && (dims[1] != 0))
  {
    if (fPrintSizeY > dims[1])
    {
      return dims[1];
    }
  }
  if (fPrintSizeY < -1)
  {
    return 0;
  }
  return fPrintSizeY;
}

void G4VulkanViewer::setExportSize(G4int X, G4int Y)
{
  fPrintSizeX = X;
  fPrintSizeY = Y;
}

/**
 If name is "" or "!", filename and extension will have the default value.
 If name is "toto.png", set the name to "toto" and the format to "png". No incremented suffix is
 added. If name is "toto", set the name to "toto" and the format to default (or current format if
 specify). If name is the same as previous, do not reset incremented suffix.
*/
bool G4VulkanViewer::setExportFilename(G4String name, G4bool inc)
{
  if (name == "!")
  {
    name = "";
  }

  if (inc)
  {
    if ((name != "") && (fExportFilename != name))
    {
      fExportFilenameIndex = 0;
    }
  }
  else
  {
    fExportFilenameIndex = -1;
  }

  if (name.size() == 0)
  {
    name = getRealPrintFilename().c_str();
  }
  else
  {
    // guess format by extention
    std::string extension = name.substr(name.find_last_of(".") + 1);
    // If there is a dot in the name the above might find rubbish, so...
    if (extension.size() >= 3 && extension.size() <= 4)
    {  // Possible extension
      if (setExportImageFormat(extension, false))
      {  // Extension found
        fExportFilename = name.substr(0, name.find_last_of("."));
      }
      else
      {  // No viable extension found
        return false;
      }
    }
    else
    {  // Assume name is already the required without-extension part
      fExportFilename = name;
    }
  }
  return true;
}

std::string G4VulkanViewer::getRealPrintFilename()
{
  std::string temp = fExportFilename;
  if (fExportFilenameIndex != -1)
  {
    temp += std::string("_");
    std::ostringstream os;
    os << std::setw(4) << std::setfill('0') << fExportFilenameIndex;
    std::string nb_str = os.str();
    temp += nb_str;
  }
  temp += "." + fExportImageFormat;
  return temp;
}

G4double G4VulkanViewer::GetSceneNearWidth()
{
  if (!fSceneHandler.GetScene())
  {
    return 0;
  }
  const G4Point3D targetPoint =
    fSceneHandler.GetScene()->GetStandardTargetPoint() + fVP.GetCurrentTargetPoint();
  G4double radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
  if (radius <= 0.) radius = 1.;
  const G4double cameraDistance = fVP.GetCameraDistance(radius);
  const VKdouble pnear = fVP.GetNearDistance(cameraDistance, radius);
  return 2 * fVP.GetFrontHalfHeight(pnear, radius);
}

VKdouble G4VulkanViewer::getSceneFarWidth()
{
  if (!fSceneHandler.GetScene())
  {
    return 0;
  }
  const G4Point3D targetPoint =
    fSceneHandler.GetScene()->GetStandardTargetPoint() + fVP.GetCurrentTargetPoint();
  G4double radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
  if (radius <= 0.) radius = 1.;
  const G4double cameraDistance = fVP.GetCameraDistance(radius);
  const VKdouble pnear = fVP.GetNearDistance(cameraDistance, radius);
  const VKdouble pfar = fVP.GetFarDistance(cameraDistance, pnear, radius);
  return 2 * fVP.GetFrontHalfHeight(pfar, radius);
}

VKdouble G4VulkanViewer::getSceneDepth()
{
  if (!fSceneHandler.GetScene())
  {
    return 0;
  }
  const G4Point3D targetPoint =
    fSceneHandler.GetScene()->GetStandardTargetPoint() + fVP.GetCurrentTargetPoint();
  G4double radius = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
  if (radius <= 0.) radius = 1.;
  const G4double cameraDistance = fVP.GetCameraDistance(radius);
  const VKdouble pnear = fVP.GetNearDistance(cameraDistance, radius);
  return fVP.GetFarDistance(cameraDistance, pnear, radius) - pnear;
}

void G4VulkanViewer::rotateScene(G4double dx, G4double dy)
{
  if (fVP.GetRotationStyle() == G4ViewParameters::freeRotation)
  {
    rotateSceneInViewDirection(dx, dy);
  }
  else
  {
    if (dx != 0)
    {
      rotateSceneThetaPhi(dx, 0);
    }
    if (dy != 0)
    {
      rotateSceneThetaPhi(0, dy);
    }
  }
}

void G4VulkanViewer::rotateSceneToggle(G4double dx, G4double dy)
{
  if (fVP.GetRotationStyle() != G4ViewParameters::freeRotation)
  {
    rotateSceneInViewDirection(dx, dy);
  }
  else
  {
    if (dx != 0)
    {
      rotateSceneThetaPhi(dx, 0);
    }
    if (dy != 0)
    {
      rotateSceneThetaPhi(0, dy);
    }
  }
}

void G4VulkanViewer::rotateSceneThetaPhi(G4double dx, G4double dy)
{
  if (!fSceneHandler.GetScene())
  {
    return;
  }

  G4Vector3D vp;
  G4Vector3D up;

  G4Vector3D xprime;
  G4Vector3D yprime;
  G4Vector3D zprime;

  G4double delta_alpha;
  G4double delta_theta;

  G4Vector3D new_vp;
  G4Vector3D new_up;

  G4double cosalpha;
  G4double sinalpha;

  G4Vector3D a1;
  G4Vector3D a2;
  G4Vector3D delta;
  G4Vector3D viewPoint;

  // phi spin stuff here

  vp = fVP.GetViewpointDirection().unit();
  up = fVP.GetUpVector().unit();

  yprime = (up.cross(vp)).unit();
  zprime = (vp.cross(yprime)).unit();

  if (fVP.GetLightsMoveWithCamera())
  {
    delta_alpha = dy * fRot_sens;
    delta_theta = -dx * fRot_sens;
  }
  else
  {
    delta_alpha = -dy * fRot_sens;
    delta_theta = dx * fRot_sens;
  }

  delta_alpha *= CLHEP::deg;
  delta_theta *= CLHEP::deg;

  new_vp = std::cos(delta_alpha) * vp + std::sin(delta_alpha) * zprime;

  // to avoid z rotation flipping
  // to allow more than 360° rotation

  if (fVP.GetLightsMoveWithCamera())
  {
    new_up = (new_vp.cross(yprime)).unit();
    if (new_vp.z() * vp.z() < 0)
    {
      new_up.set(new_up.x(), -new_up.y(), new_up.z());
    }
  }
  else
  {
    new_up = up;
    if (new_vp.z() * vp.z() < 0)
    {
      new_up.set(new_up.x(), -new_up.y(), new_up.z());
    }
  }
  fVP.SetUpVector(new_up);
  ////////////////
  // Rotates by fixed azimuthal angle delta_theta.

  cosalpha = new_up.dot(new_vp.unit());
  sinalpha = std::sqrt(1. - std::pow(cosalpha, 2));
  yprime = (new_up.cross(new_vp.unit())).unit();
  xprime = yprime.cross(new_up);
  // Projection of vp on plane perpendicular to up...
  a1 = sinalpha * xprime;
  // Required new projection...
  a2 = sinalpha * (std::cos(delta_theta) * xprime + std::sin(delta_theta) * yprime);
  // Required Increment vector...
  delta = a2 - a1;
  // So new viewpoint is...
  viewPoint = new_vp.unit() + delta;

  fVP.SetViewAndLights(viewPoint);
}

void G4VulkanViewer::rotateSceneInViewDirection(G4double dx, G4double dy)
{
  if (!fSceneHandler.GetScene())
  {
    return;
  }

  G4Vector3D vp;
  G4Vector3D up;

  G4Vector3D xprime;
  G4Vector3D yprime;
  G4Vector3D zprime;

  G4Vector3D new_vp;
  G4Vector3D new_up;

  G4Vector3D a1;
  G4Vector3D a2;
  G4Vector3D delta;
  G4Vector3D viewPoint;

  dx = dx / 100;
  dy = dy / 100;

  // phi spin stuff here

  vp = fVP.GetViewpointDirection().unit();
  up = fVP.GetUpVector().unit();

  G4Vector3D zPrimeVector =
    G4Vector3D(up.y() * vp.z() - up.z() * vp.y(), up.z() * vp.x() - up.x() * vp.z(),
               up.x() * vp.y() - up.y() * vp.x());

  viewPoint = vp / fRot_sens + (zPrimeVector * dx - up * dy);
  new_up = G4Vector3D(viewPoint.y() * zPrimeVector.z() - viewPoint.z() * zPrimeVector.y(),
                      viewPoint.z() * zPrimeVector.x() - viewPoint.x() * zPrimeVector.z(),
                      viewPoint.x() * zPrimeVector.y() - viewPoint.y() * zPrimeVector.x());

  G4Vector3D new_upUnit = new_up.unit();

  fVP.SetUpVector(new_upUnit);
  fVP.SetViewAndLights(viewPoint);
}

void G4VulkanViewer::addExportImageFormat(std::string format)
{
  fExportImageFormatVector.push_back(format);
}

bool G4VulkanViewer::setExportImageFormat(std::string format, bool quiet)
{
  bool found = false;
  std::string list;
  for (unsigned int a = 0; a < fExportImageFormatVector.size(); a++)
  {
    list += fExportImageFormatVector.at(a) + " ";

    if (fExportImageFormatVector.at(a) == format)
    {
      if (!quiet)
      {
        G4cout << " Changing export format to \"" << format << "\"" << G4endl;
      }
      if (format != fExportImageFormat)
      {
        fExportFilenameIndex = 0;
        fExportImageFormat = format;
      }
      return true;
    }
  }
  if (!found)
  {
    if (format.size() == 0)
    {
      G4cout << " Current formats availables are : " << list << G4endl;
    }
    else
    {
      G4cerr << " Format \"" << format
             << "\" is not available for the selected viewer. Current formats availables are : "
             << list << G4endl;
    }
  }
  return false;
}

// From MESA implementation :
// http://www.techques.com/question/1-8660454/gluPickMatrix-code-from-Mesa

void G4VulkanViewer::g4VkPickMatrix(VKdouble x, VKdouble y, VKdouble width, VKdouble height,
                                     VKint viewport[4])
{
  VKdouble mat[16];
  VKdouble sx, sy;
  VKdouble tx, ty;

  sx = viewport[2] / width;
  sy = viewport[3] / height;
  tx = (viewport[2] + 2.0 * (viewport[0] - x)) / width;
  ty = (viewport[3] + 2.0 * (viewport[1] - y)) / height;

#define M(row, col) mat[col * 4 + row]
  M(0, 0) = sx;
  M(0, 1) = 0.0;
  M(0, 2) = 0.0;
  M(0, 3) = tx;
  M(1, 0) = 0.0;
  M(1, 1) = sy;
  M(1, 2) = 0.0;
  M(1, 3) = ty;
  M(2, 0) = 0.0;
  M(2, 1) = 0.0;
  M(2, 2) = 1.0;
  M(2, 3) = 0.0;
  M(3, 0) = 0.0;
  M(3, 1) = 0.0;
  M(3, 2) = 0.0;
  M(3, 3) = 1.0;
#undef M

  vkMultMatrixd(mat);
}

// From MESA implementation :
// https://github.com/jlamarche/iOS-OpenGLES-Stuff/blob/master/Wavefront%20OBJ%20Loader/Classes/gluLookAt.m
// or
// http://www.daniweb.com/software-development/game-development/threads/308901/lookat-matrix-source-code

void G4VulkanViewer::g4VkLookAt(VKdouble eyex, VKdouble eyey, VKdouble eyez, VKdouble centerx,
                                 VKdouble centery, VKdouble centerz, VKdouble upx, VKdouble upy,
                                 VKdouble upz)
{
  VKdouble mat[16];
  VKdouble x[3], y[3], z[3];
  VKdouble mag;

  /* Make rotation matrix */

  /* Z vector */
  z[0] = eyex - centerx;
  z[1] = eyey - centery;
  z[2] = eyez - centerz;
  mag = std::sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);
  if (mag)
  { /* mpichler, 19950515 */
    z[0] /= mag;
    z[1] /= mag;
    z[2] /= mag;
  }

  /* Y vector */
  y[0] = upx;
  y[1] = upy;
  y[2] = upz;

  /* X vector = Y cross Z */
  x[0] = y[1] * z[2] - y[2] * z[1];
  x[1] = -y[0] * z[2] + y[2] * z[0];
  x[2] = y[0] * z[1] - y[1] * z[0];

  /* Recompute Y = Z cross X */
  y[0] = z[1] * x[2] - z[2] * x[1];
  y[1] = -z[0] * x[2] + z[2] * x[0];
  y[2] = z[0] * x[1] - z[1] * x[0];

  /* mpichler, 19950515 */
  /* cross product gives area of parallelogram, which is < 1.0 for
   * non-perpendicular unit-length vectors; so normalize x, y here
   */

  mag = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  if (mag)
  {
    x[0] /= mag;
    x[1] /= mag;
    x[2] /= mag;
  }

  mag = std::sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
  if (mag)
  {
    y[0] /= mag;
    y[1] /= mag;
    y[2] /= mag;
  }

#define M(row, col) mat[col * 4 + row]
  M(0, 0) = x[0];
  M(0, 1) = x[1];
  M(0, 2) = x[2];
  M(0, 3) = 0.0;
  M(1, 0) = y[0];
  M(1, 1) = y[1];
  M(1, 2) = y[2];
  M(1, 3) = 0.0;
  M(2, 0) = z[0];
  M(2, 1) = z[1];
  M(2, 2) = z[2];
  M(2, 3) = 0.0;
  M(3, 0) = 0.0;
  M(3, 1) = 0.0;
  M(3, 2) = 0.0;
  M(3, 3) = 1.0;
#undef M
  vkMultMatrixd(mat);

  /* Translate Eye to Origin */
  vkTranslated(-eyex, -eyey, -eyez);
}

void G4VulkanViewer::g4VkOrtho(VKdouble left, VKdouble right, VKdouble bottom, VKdouble top,
                               VKdouble zNear, VKdouble zFar)
{
  //  vkOrtho (left, right, bottom, top, near, far);

  VKdouble a = 2.0 / (right - left);
  VKdouble b = 2.0 / (top - bottom);
  VKdouble c = -2.0 / (zFar - zNear);

  VKdouble tx = -(right + left) / (right - left);
  VKdouble ty = -(top + bottom) / (top - bottom);
  VKdouble tz = -(zFar + zNear) / (zFar - zNear);

  VKdouble ortho[16] = {a, 0, 0, 0, 0, b, 0, 0, 0, 0, c, 0, tx, ty, tz, 1};
  vkMultMatrixd(ortho);
}

void G4VulkanViewer::g4VkFrustum(VKdouble left, VKdouble right, VKdouble bottom, VKdouble top,
                                 VKdouble zNear, VKdouble zFar)
{
  //  vkFrustum (left, right, bottom, top, near, far);

  VKdouble deltaX = right - left;
  VKdouble deltaY = top - bottom;
  VKdouble deltaZ = zFar - zNear;

  VKdouble a = 2.0f * zNear / deltaX;
  VKdouble b = 2.0f * zNear / deltaY;
  VKdouble c = (right + left) / deltaX;
  VKdouble d = (top + bottom) / deltaY;
  VKdouble e = -(zFar + zNear) / (zFar - zNear);
  VKdouble f = -2.0f * zFar * zNear / deltaZ;

  VKdouble proj[16] = {a, 0, 0, 0, 0, b, 0, 0, c, d, e, -1.0f, 0, 0, f, 0};

  vkMultMatrixd(proj);
}

G4String G4VulkanViewerPickMap::print()
{
  std::ostringstream txt;
  for (unsigned int a = 0; a < fAttributes.size(); a++)
  {
    txt << fAttributes[a];
    if (a < fAttributes.size() - 1) txt << "\n";
  }
  return txt.str();
}
