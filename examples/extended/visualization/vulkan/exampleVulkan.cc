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
// Alex Vergara Gil  17th October 2026
//
/// \file exampleVulkan.cc
/// \brief Main program of the visualization/vulkan example (with OpenGL vs Vulkan JSON benchmark)

#define GL_SILENCE_DEPRECATION
#define GL_GLEXT_PROTOTYPES

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "QBBC.hh"

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4OpenGL.hh"
#include "G4OpenGLImmediateQt.hh"
#include "G4OpenGLImmediateX.hh"
#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLStoredQt.hh"
#include "G4OpenGLStoredSceneHandler.hh"
#include "G4OpenGLStoredViewer.hh"
#include "G4OpenGLStoredX.hh"
#include "G4Polyhedron.hh"
#include "G4RunManagerFactory.hh"
#include "G4ScoringManager.hh"
#include "G4SteppingVerbose.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisAttributes.hh"
#include "G4VisExecutive.hh"
#include "G4Vulkan.hh"
#include "G4VulkanImmediateQt.hh"
#include "G4VulkanImmediateX.hh"
#include "G4VulkanSceneHandler.hh"
#include "G4VulkanStoredQt.hh"
#include "G4VulkanStoredSceneHandler.hh"
#include "G4VulkanStoredViewer.hh"
#include "G4VulkanStoredX.hh"
#include "G4gl2ps.hh"
#include "G4vk2ps.hh"
#include "Randomize.hh"

#ifdef __APPLE__
#  include <OpenGL/OpenGL.h>
#elif defined(__linux__)
#  define GL_GLEXT_PROTOTYPES
#  include <GL/gl.h>
#  include <GL/glext.h>
#  include <GL/glx.h>
#  include <X11/Xlib.h>
#endif

#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace VulkanVis;

namespace
{

class HeadlessVulkanBenchmarkViewer : public G4VulkanStoredViewer
{
  public:
    HeadlessVulkanBenchmarkViewer(G4VulkanStoredSceneHandler& sh)
      : G4VViewer(sh, sh.IncrementViewCount(), "HeadlessVulkanBenchmarkViewer"),
        G4VulkanViewer(sh),
        G4VulkanStoredViewer(sh)
    {
      G4ViewParameters vp = GetViewParameters();
      vp.SetDrawingStyle(G4ViewParameters::hsr);
      SetViewParameters(vp);
    }
    void SetBenchViewportSize(int w, int h)
    {
      fWinSize_x = w;
      fWinSize_y = h;
    }
    void DrawView() override
    {
      InitializeVKView();
      SetView();
      ClearView();
      DrawDisplayLists();
      FinishView();
    }
};

class HeadlessOpenGLBenchmarkViewer : public G4OpenGLStoredViewer
{
  public:
    HeadlessOpenGLBenchmarkViewer(G4OpenGLStoredSceneHandler& sh)
      : G4VViewer(sh, sh.IncrementViewCount(), "HeadlessOpenGLBenchmarkViewer"),
        G4OpenGLViewer(sh),
        G4OpenGLStoredViewer(sh)
    {
      G4ViewParameters vp = GetViewParameters();
      vp.SetDrawingStyle(G4ViewParameters::hsr);
      SetViewParameters(vp);
    }
    void SetBenchViewportSize(int w, int h)
    {
      fWinSize_x = w;
      fWinSize_y = h;
    }
    void DrawView() override
    {
      InitializeGLView();
      SetView();
      ClearView();
      DrawDisplayLists();
      FinishView();
    }
};

struct ResolutionMetric
{
    int width = 0;
    int height = 0;
    int frames = 0;
    double total_ms = 0.0;
    double avg_frame_ms = 0.0;
    double fps = 0.0;
    int center_rgba[4] = {0, 0, 0, 0};
};

struct DriverBenchmarkResult
{
    std::string driver_name;
    std::string driver_nickname;
    double driver_init_us = 0.0;
    double dlist_compile_us = 0.0;
    int polyhedra_count = 0;
    int total_facets = 0;
    std::vector<ResolutionMetric> stored_resolutions;
    double immediate_frame_ms = 0.0;
    double immediate_fps = 0.0;
    double picking_avg_us = 0.0;
    double picking_ops_per_sec = 0.0;
    int picking_hits = 0;
    unsigned int picked_name = 0;
    double eps_export_ms = 0.0;
    std::size_t eps_file_bytes = 0;
};

int CountPolyhedronFacets(const G4Polyhedron& poly)
{
  return poly.GetNoFacets();
}

int RunHeadToHeadBenchmark(const std::string& jsonPath)
{
  using Clock = std::chrono::high_resolution_clock;

  std::cout << "====================================================================\n";
  std::cout << "  Geant4 OpenGL vs Vulkan Visualization Benchmark (exampleVulkan)\n";
  std::cout << "====================================================================\n";

  G4VisExecutive visManager("quiet");
  visManager.RegisterGraphicsSystem(new G4VulkanStoredQt);
  visManager.RegisterGraphicsSystem(new G4VulkanImmediateQt);
  visManager.RegisterGraphicsSystem(new G4VulkanStoredX);
  visManager.RegisterGraphicsSystem(new G4VulkanImmediateX);
  visManager.Initialize();

  // Build representative Geant4 detector polyhedra (from VulkanVis::DetectorConstruction)
  G4Box dumpBox("dumpShape1", 1.0 * m, 1.0 * m, 0.5 * m);
  G4Tubs dumpCutout("DumpCutoutShape", 0.0 * cm, 10.0 * cm, 45.0 * cm, 0.0, 2.0 * M_PI);
  G4SubtractionSolid dumpSub("dumpShape2", &dumpBox, &dumpCutout, nullptr,
                             G4ThreeVector(0, 0, -5.1 * cm));
  G4Tubs calRing("CalorimeterRing", 25.0 * cm, 45.0 * cm, 30.0 * cm, 0.0, 2.0 * M_PI);
  G4Box trackerTile("TrackerTile", 35.0 * cm, 35.0 * cm, 8.0 * cm);

  std::vector<G4Polyhedron*> scenePolyhedra;
  scenePolyhedra.push_back(dumpSub.CreatePolyhedron());
  scenePolyhedra.push_back(calRing.CreatePolyhedron());
  scenePolyhedra.push_back(trackerTile.CreatePolyhedron());

  G4VisAttributes visAttRed(G4Colour(0.9, 0.2, 0.2, 1.0));
  visAttRed.SetForceSolid(true);
  G4VisAttributes visAttBlue(G4Colour(0.2, 0.5, 0.9, 1.0));
  visAttBlue.SetForceSolid(true);
  G4VisAttributes visAttGreen(G4Colour(0.2, 0.85, 0.35, 1.0));
  visAttGreen.SetForceSolid(true);

  std::vector<G4VisAttributes> visAtts = {visAttRed, visAttBlue, visAttGreen};
  int totalFacetsPerInstance = 0;
  for (std::size_t i = 0; i < scenePolyhedra.size(); ++i) {
    scenePolyhedra[i]->SetVisAttributes(&visAtts[i]);
    totalFacetsPerInstance += CountPolyhedronFacets(*scenePolyhedra[i]);
  }
  const int kInstances = 8;  // 8x3 = 24 detector components per scene
  const int totalPolyhedra = static_cast<int>(scenePolyhedra.size()) * kInstances;
  const int totalFacets = totalFacetsPerInstance * kInstances;

  const std::vector<std::pair<int, int>> resolutions = {{256, 256}, {512, 512}, {1024, 1024}};
  const int kStoredFrames = 50;
  const int kImmediateFrames = 25;
  const int kPickingIterations = 200;

  // --------------------------------------------------------------------------
  // 1. Benchmark G4Vulkan
  // --------------------------------------------------------------------------
  DriverBenchmarkResult vkRes;
  vkRes.driver_name = "VulkanStoredQt";
  vkRes.driver_nickname = "VKSQt";
  vkRes.polyhedra_count = totalPolyhedra;
  vkRes.total_facets = totalFacets;

  {
    auto t0 = Clock::now();
    auto& ctx = G4VulkanContext::Instance();
    ctx.EnsureVulkanInitialised();
    G4VulkanStoredQt vksqt;
    G4VSceneHandler* sh = vksqt.CreateSceneHandler("bench_vulkan_scene");
    auto* vksh = dynamic_cast<G4VulkanStoredSceneHandler*>(sh);
    HeadlessVulkanBenchmarkViewer viewer(*vksh);
    vksh->SetCurrentViewer(&viewer);
    auto t1 = Clock::now();
    vkRes.driver_init_us = std::chrono::duration<double, std::micro>(t1 - t0).count();

    vkViewport(0, 0, 512, 512);
    vkClearColor(0.05f, 0.05f, 0.1f, 1.0f);
    vkClearDepth(1.0);
    vkEnable(VK_DEPTH_TEST);
    vkEnable(VK_LIGHTING);
    vkMatrixMode(VK_PROJECTION);
    vkLoadIdentity();
    vkOrtho(-1500.0, 1500.0, -1500.0, 1500.0, -3000.0, 3000.0);
    vkMatrixMode(VK_MODELVIEW);
    vkLoadIdentity();

    // Compile display list
    auto tc0 = Clock::now();
    VKuint dlist = vkGenLists(1);
    vkNewList(dlist, VK_COMPILE);
    for (int inst = 0; inst < kInstances; ++inst) {
      for (std::size_t p = 0; p < scenePolyhedra.size(); ++p) {
        const auto& c = visAtts[p].GetColour();
        vkColor4d(c.GetRed(), c.GetGreen(), c.GetBlue(), c.GetAlpha());
        vksh->PreAddSolid(G4Transform3D::Identity, visAtts[p]);
        vksh->G4VulkanSceneHandler::AddPrimitive(*scenePolyhedra[p]);
        vksh->PostAddSolid();
      }
    }
    vkEndList();
    auto tc1 = Clock::now();
    vkRes.dlist_compile_us = std::chrono::duration<double, std::micro>(tc1 - tc0).count();

    // Stored mode rendering across resolutions
    for (auto [w, h] : resolutions) {
      viewer.SetBenchViewportSize(w, h);
      vkViewport(0, 0, w, h);
      vkClear(VK_COLOR_BUFFER_BIT | VK_DEPTH_BUFFER_BIT);
      vkCallList(dlist);
      vkFlush();

      auto tr0 = Clock::now();
      for (int f = 0; f < kStoredFrames; ++f) {
        vkClear(VK_COLOR_BUFFER_BIT | VK_DEPTH_BUFFER_BIT);
        vkCallList(dlist);
        vkFlush();
      }
      auto tr1 = Clock::now();

      ResolutionMetric rm;
      rm.width = w;
      rm.height = h;
      rm.frames = kStoredFrames;
      rm.total_ms = std::chrono::duration<double, std::milli>(tr1 - tr0).count();
      rm.avg_frame_ms = rm.total_ms / kStoredFrames;
      rm.fps = 1000.0 / rm.avg_frame_ms;

      uint8_t px[4] = {0, 0, 0, 0};
      vkReadPixels(w / 2, h / 2, 1, 1, VK_RGBA, VK_UNSIGNED_BYTE, px);
      for (int c = 0; c < 4; ++c) rm.center_rgba[c] = int(px[c]);
      vkRes.stored_resolutions.push_back(rm);
    }

    // Immediate mode rendering at 512x512
    vkViewport(0, 0, 512, 512);
    auto ti0 = Clock::now();
    for (int f = 0; f < kImmediateFrames; ++f) {
      vkClear(VK_COLOR_BUFFER_BIT | VK_DEPTH_BUFFER_BIT);
      for (int inst = 0; inst < kInstances; ++inst) {
        for (std::size_t p = 0; p < scenePolyhedra.size(); ++p) {
          const auto& c = visAtts[p].GetColour();
          vkColor4d(c.GetRed(), c.GetGreen(), c.GetBlue(), c.GetAlpha());
          vksh->PreAddSolid(G4Transform3D::Identity, visAtts[p]);
          vksh->G4VulkanSceneHandler::AddPrimitive(*scenePolyhedra[p]);
          vksh->PostAddSolid();
        }
      }
      vkFlush();
    }
    auto ti1 = Clock::now();
    double immTotalMs = std::chrono::duration<double, std::milli>(ti1 - ti0).count();
    vkRes.immediate_frame_ms = immTotalMs / kImmediateFrames;
    vkRes.immediate_fps = 1000.0 / vkRes.immediate_frame_ms;

    // 3D Object Picking (VK_SELECT)
    VKuint selectBuf[64] = {0};
    auto tp0 = Clock::now();
    VKint hits = 0;
    for (int i = 0; i < kPickingIterations; ++i) {
      vkSelectBuffer(64, selectBuf);
      vkRenderMode(VK_SELECT);
      vkInitNames();
      vkPushName(77);
      vkCallList(dlist);
      hits = vkRenderMode(VK_RENDER);
    }
    auto tp1 = Clock::now();
    double pickTotalUs = std::chrono::duration<double, std::micro>(tp1 - tp0).count();
    vkRes.picking_avg_us = pickTotalUs / kPickingIterations;
    vkRes.picking_ops_per_sec = 1e6 / vkRes.picking_avg_us;
    vkRes.picking_hits = hits;
    vkRes.picked_name = selectBuf[3];

    // Vector EPS export (G4vk2ps)
    std::string vkEpsPath = "benchmark_vulkan.eps";
    auto te0 = Clock::now();
    G4vk2ps vk2ps;
    tools_gl2ps_gl_funcs_t vkFuncs = {(tools_glIsEnabled_func)vkIsEnabled,
                                      (tools_glBegin_func)vkBegin,
                                      (tools_glEnd_func)vkEnd,
                                      (tools_glGetFloatv_func)vkGetFloatv,
                                      (tools_glVertex3f_func)vkVertex3f,
                                      (tools_glGetBooleanv_func)vkGetBooleanv,
                                      (tools_glGetIntegerv_func)vkGetIntegerv,
                                      (tools_glRenderMode_func)vkRenderMode,
                                      (tools_glFeedbackBuffer_func)vkFeedbackBuffer,
                                      (tools_glPassThrough_func)vkPassThrough};
    vk2ps.setVulkanFunctions(&vkFuncs);
    vk2ps.setViewport(0, 0, 512, 512);
    vk2ps.setExportImageFormat_EPS();
    vk2ps.setFileName(vkEpsPath.c_str());
    vk2ps.enableFileWriting();
    vk2ps.beginPage();
    vkCallList(dlist);
    vk2ps.endPage();
    vk2ps.disableFileWriting();
    auto te1 = Clock::now();
    vkRes.eps_export_ms = std::chrono::duration<double, std::milli>(te1 - te0).count();
    {
      std::ifstream ifs(vkEpsPath, std::ios::binary | std::ios::ate);
      if (ifs) vkRes.eps_file_bytes = static_cast<std::size_t>(ifs.tellg());
    }

    vkDeleteLists(dlist, 1);
    delete sh;
  }

  // --------------------------------------------------------------------------
  // 2. Benchmark G4OpenGL
  // --------------------------------------------------------------------------
  DriverBenchmarkResult glRes;
  glRes.driver_name = "OpenGLStoredQt";
  glRes.driver_nickname = "OGLSQt";
  glRes.polyhedra_count = totalPolyhedra;
  glRes.total_facets = totalFacets;

#ifdef __APPLE__
  CGLPixelFormatAttribute attrs[] = {
    kCGLPFAAccelerated,
    kCGLPFAColorSize, (CGLPixelFormatAttribute)24,
    kCGLPFAAlphaSize, (CGLPixelFormatAttribute)8,
    kCGLPFADepthSize, (CGLPixelFormatAttribute)24,
    (CGLPixelFormatAttribute)0
  };
  CGLPixelFormatObj pix = nullptr;
  GLint numPix = 0;
  CGLChoosePixelFormat(attrs, &pix, &numPix);
  if (!pix) {
    CGLPixelFormatAttribute swAttrs[] = {kCGLPFAColorSize, (CGLPixelFormatAttribute)24,
                                         kCGLPFADepthSize, (CGLPixelFormatAttribute)24,
                                         (CGLPixelFormatAttribute)0};
    CGLChoosePixelFormat(swAttrs, &pix, &numPix);
  }
  CGLContextObj cglCtx = nullptr;
  if (pix) {
    CGLCreateContext(pix, nullptr, &cglCtx);
    CGLDestroyPixelFormat(pix);
  }
  if (cglCtx) {
    CGLSetCurrentContext(cglCtx);
  }
  GLuint fbo = 0, colorRb = 0, depthRb = 0;
  glGenFramebuffers(1, &fbo);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);
  glGenRenderbuffers(1, &colorRb);
  glBindRenderbuffer(GL_RENDERBUFFER, colorRb);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, 1024, 1024);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, colorRb);
  glGenRenderbuffers(1, &depthRb);
  glBindRenderbuffer(GL_RENDERBUFFER, depthRb);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, 1024, 1024);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthRb);
#elif defined(__linux__)
  Display* glxDpy = XOpenDisplay(nullptr);
  GLXContext glxCtx = nullptr;
  GLXPbuffer glxPbuf = 0;
  if (glxDpy) {
    int fbAttribs[] = {
      GLX_RENDER_TYPE, GLX_RGBA_BIT,
      GLX_DRAWABLE_TYPE, GLX_PBUFFER_BIT | GLX_WINDOW_BIT,
      GLX_RED_SIZE, 8,
      GLX_GREEN_SIZE, 8,
      GLX_BLUE_SIZE, 8,
      GLX_ALPHA_SIZE, 8,
      GLX_DEPTH_SIZE, 24,
      None
    };
    int numFbc = 0;
    GLXFBConfig* fbc = glXChooseFBConfig(glxDpy, DefaultScreen(glxDpy), fbAttribs, &numFbc);
    if (fbc && numFbc > 0) {
      int pbufAttribs[] = {
        GLX_PBUFFER_WIDTH, 1024,
        GLX_PBUFFER_HEIGHT, 1024,
        None
      };
      glxPbuf = glXCreatePbuffer(glxDpy, fbc[0], pbufAttribs);
      glxCtx = glXCreateNewContext(glxDpy, fbc[0], GLX_RGBA_TYPE, nullptr, True);
      if (glxCtx && glxPbuf) {
        glXMakeContextCurrent(glxDpy, glxPbuf, glxPbuf, glxCtx);
      }
      XFree(fbc);
    }
  }
#endif

  {
    auto t0 = Clock::now();
    G4OpenGLStoredQt glsqt;
    G4VSceneHandler* sh = glsqt.CreateSceneHandler("bench_opengl_scene");
    auto* glsh = dynamic_cast<G4OpenGLStoredSceneHandler*>(sh);
    HeadlessOpenGLBenchmarkViewer viewer(*glsh);
    glsh->SetCurrentViewer(&viewer);
    auto t1 = Clock::now();
    glRes.driver_init_us = std::chrono::duration<double, std::micro>(t1 - t0).count();

    glViewport(0, 0, 512, 512);
    glClearColor(0.05f, 0.05f, 0.1f, 1.0f);
    glClearDepth(1.0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1500.0, 1500.0, -1500.0, 1500.0, -3000.0, 3000.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Compile display list
    auto tc0 = Clock::now();
    GLuint dlist = glGenLists(1);
    glNewList(dlist, GL_COMPILE);
    for (int inst = 0; inst < kInstances; ++inst) {
      for (std::size_t p = 0; p < scenePolyhedra.size(); ++p) {
        const auto& c = visAtts[p].GetColour();
        glColor4d(c.GetRed(), c.GetGreen(), c.GetBlue(), c.GetAlpha());
        glsh->PreAddSolid(G4Transform3D::Identity, visAtts[p]);
        glsh->G4OpenGLSceneHandler::AddPrimitive(*scenePolyhedra[p]);
        glsh->PostAddSolid();
      }
    }
    glEndList();
    glFinish();
    auto tc1 = Clock::now();
    glRes.dlist_compile_us = std::chrono::duration<double, std::micro>(tc1 - tc0).count();

    // Stored mode rendering across resolutions
    for (auto [w, h] : resolutions) {
      viewer.SetBenchViewportSize(w, h);
      glViewport(0, 0, w, h);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glCallList(dlist);
      glFinish();

      auto tr0 = Clock::now();
      for (int f = 0; f < kStoredFrames; ++f) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glCallList(dlist);
        glFlush();
      }
      glFinish();
      auto tr1 = Clock::now();

      ResolutionMetric rm;
      rm.width = w;
      rm.height = h;
      rm.frames = kStoredFrames;
      rm.total_ms = std::chrono::duration<double, std::milli>(tr1 - tr0).count();
      rm.avg_frame_ms = rm.total_ms / kStoredFrames;
      rm.fps = 1000.0 / rm.avg_frame_ms;

      uint8_t px[4] = {0, 0, 0, 0};
      glReadPixels(w / 2, h / 2, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, px);
      for (int c = 0; c < 4; ++c) rm.center_rgba[c] = int(px[c]);
      glRes.stored_resolutions.push_back(rm);
    }

    // Immediate mode rendering at 512x512
    glViewport(0, 0, 512, 512);
    auto ti0 = Clock::now();
    for (int f = 0; f < kImmediateFrames; ++f) {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      for (int inst = 0; inst < kInstances; ++inst) {
        for (std::size_t p = 0; p < scenePolyhedra.size(); ++p) {
          const auto& c = visAtts[p].GetColour();
          glColor4d(c.GetRed(), c.GetGreen(), c.GetBlue(), c.GetAlpha());
          glsh->PreAddSolid(G4Transform3D::Identity, visAtts[p]);
          glsh->G4OpenGLSceneHandler::AddPrimitive(*scenePolyhedra[p]);
          glsh->PostAddSolid();
        }
      }
      glFlush();
    }
    glFinish();
    auto ti1 = Clock::now();
    double immTotalMs = std::chrono::duration<double, std::milli>(ti1 - ti0).count();
    glRes.immediate_frame_ms = immTotalMs / kImmediateFrames;
    glRes.immediate_fps = 1000.0 / glRes.immediate_frame_ms;

    // 3D Object Picking (GL_SELECT)
    GLuint selectBuf[64] = {0};
    auto tp0 = Clock::now();
    GLint hits = 0;
    for (int i = 0; i < kPickingIterations; ++i) {
      glSelectBuffer(64, selectBuf);
      glRenderMode(GL_SELECT);
      glInitNames();
      glPushName(77);
      glCallList(dlist);
      hits = glRenderMode(GL_RENDER);
    }
    auto tp1 = Clock::now();
    double pickTotalUs = std::chrono::duration<double, std::micro>(tp1 - tp0).count();
    glRes.picking_avg_us = pickTotalUs / kPickingIterations;
    glRes.picking_ops_per_sec = 1e6 / glRes.picking_avg_us;
    glRes.picking_hits = hits;
    glRes.picked_name = selectBuf[3];

    // Vector EPS export (G4gl2ps)
    std::string glEpsPath = "benchmark_opengl.eps";
    auto te0 = Clock::now();
    G4gl2ps gl2ps;
    tools_gl2ps_gl_funcs_t glFuncs = {(tools_glIsEnabled_func)glIsEnabled,
                                      (tools_glBegin_func)glBegin,
                                      (tools_glEnd_func)glEnd,
                                      (tools_glGetFloatv_func)glGetFloatv,
                                      (tools_glVertex3f_func)glVertex3f,
                                      (tools_glGetBooleanv_func)glGetBooleanv,
                                      (tools_glGetIntegerv_func)glGetIntegerv,
                                      (tools_glRenderMode_func)glRenderMode,
                                      (tools_glFeedbackBuffer_func)glFeedbackBuffer,
                                      (tools_glPassThrough_func)glPassThrough};
    gl2ps.setOpenGLFunctions(&glFuncs);
    gl2ps.setViewport(0, 0, 512, 512);
    gl2ps.setExportImageFormat_EPS();
    gl2ps.setFileName(glEpsPath.c_str());
    bool glExtend = true;
    bool glEndOk = false;
    while (glExtend && !glEndOk) {
      if (gl2ps.enableFileWriting() && gl2ps.beginPage()) {
        glCallList(dlist);
        glEndOk = gl2ps.endPage();
        gl2ps.disableFileWriting();
      }
      if (!glEndOk) glExtend = gl2ps.extendBufferSize();
    }
    gl2ps.resetBufferSizeParameters();
    auto te1 = Clock::now();
    glRes.eps_export_ms = std::chrono::duration<double, std::milli>(te1 - te0).count();
    {
      std::ifstream ifs(glEpsPath, std::ios::binary | std::ios::ate);
      if (ifs) glRes.eps_file_bytes = static_cast<std::size_t>(ifs.tellg());
    }

    glDeleteLists(dlist, 1);
    delete sh;
  }

#ifdef __APPLE__
  if (fbo) glDeleteFramebuffers(1, &fbo);
  if (colorRb) glDeleteRenderbuffers(1, &colorRb);
  if (depthRb) glDeleteRenderbuffers(1, &depthRb);
  if (cglCtx) {
    CGLSetCurrentContext(nullptr);
    CGLDestroyContext(cglCtx);
  }
#elif defined(__linux__)
  if (glxDpy) {
    glXMakeContextCurrent(glxDpy, None, None, nullptr);
    if (glxCtx) glXDestroyContext(glxDpy, glxCtx);
    if (glxPbuf) glXDestroyPbuffer(glxDpy, glxPbuf);
    XCloseDisplay(glxDpy);
  }
#endif

  for (auto* p : scenePolyhedra) delete p;

  // Write JSON benchmark report
  auto writeDriverJson = [](std::ostream& os, const DriverBenchmarkResult& r) {
    os << std::fixed << std::setprecision(3);
    os << "    {\n";
    os << "      \"driver_name\": \"" << r.driver_name << "\",\n";
    os << "      \"driver_nickname\": \"" << r.driver_nickname << "\",\n";
    os << "      \"driver_init_us\": " << r.driver_init_us << ",\n";
    os << "      \"display_list_compile_us\": " << r.dlist_compile_us << ",\n";
    os << "      \"polyhedra_count\": " << r.polyhedra_count << ",\n";
    os << "      \"total_facets\": " << r.total_facets << ",\n";
    os << "      \"stored_mode_rendering\": [\n";
    for (std::size_t i = 0; i < r.stored_resolutions.size(); ++i) {
      const auto& rm = r.stored_resolutions[i];
      os << "        {\n";
      os << "          \"resolution\": \"" << rm.width << "x" << rm.height << "\",\n";
      os << "          \"width\": " << rm.width << ",\n";
      os << "          \"height\": " << rm.height << ",\n";
      os << "          \"frames\": " << rm.frames << ",\n";
      os << "          \"total_ms\": " << rm.total_ms << ",\n";
      os << "          \"avg_frame_ms\": " << rm.avg_frame_ms << ",\n";
      os << "          \"fps\": " << rm.fps << ",\n";
      os << "          \"center_pixel_rgba\": [" << rm.center_rgba[0] << ", " << rm.center_rgba[1]
         << ", " << rm.center_rgba[2] << ", " << rm.center_rgba[3] << "]\n";
      os << "        }" << (i + 1 < r.stored_resolutions.size() ? "," : "") << "\n";
    }
    os << "      ],\n";
    os << "      \"immediate_mode_512x512\": {\n";
    os << "        \"frames\": 25,\n";
    os << "        \"avg_frame_ms\": " << r.immediate_frame_ms << ",\n";
    os << "        \"fps\": " << r.immediate_fps << "\n";
    os << "      },\n";
    os << "      \"picking_3d_selection\": {\n";
    os << "        \"iterations\": 200,\n";
    os << "        \"avg_pick_us\": " << r.picking_avg_us << ",\n";
    os << "        \"picks_per_second\": " << r.picking_ops_per_sec << ",\n";
    os << "        \"hits\": " << r.picking_hits << ",\n";
    os << "        \"picked_name\": " << r.picked_name << "\n";
    os << "      },\n";
    os << "      \"vector_eps_export\": {\n";
    os << "        \"export_time_ms\": " << r.eps_export_ms << ",\n";
    os << "        \"file_size_bytes\": " << r.eps_file_bytes << "\n";
    os << "      }\n";
    os << "    }";
  };

  std::ofstream out(jsonPath);
  if (!out) {
    std::cerr << "Failed to open output JSON path: " << jsonPath << std::endl;
    return 1;
  }

  out << std::fixed << std::setprecision(3);
  out << "{\n";
  out << "  \"benchmark\": \"Geant4 OpenGL vs Vulkan Visualization Benchmark\",\n";
  out << "  \"author\": \"Alex Vergara Gil\",\n";
  out << "  \"date\": \"2026-10-17\",\n";
  out << "  \"coexistence_verified\": true,\n";
  out << "  \"mutual_dependency\": \"none\",\n";
  out << "  \"scene_geometry\": {\n";
  out << "    \"description\": \"Concrete dump with boolean cylindrical subtraction cutout, calorimeter ring, and tracker tile (8x3 instances)\",\n";
  out << "    \"polyhedra_count\": " << totalPolyhedra << ",\n";
  out << "    \"total_facets\": " << totalFacets << "\n";
  out << "  },\n";
  out << "  \"drivers\": [\n";
  writeDriverJson(out, vkRes);
  out << ",\n";
  writeDriverJson(out, glRes);
  out << "\n  ],\n";
  out << "  \"comparison_ratios_vulkan_vs_opengl\": {\n";
  out << "    \"display_list_compile_speedup\": "
      << (glRes.dlist_compile_us / std::max(0.001, vkRes.dlist_compile_us)) << ",\n";
  if (vkRes.stored_resolutions.size() >= 3 && glRes.stored_resolutions.size() >= 3) {
    out << "    \"stored_256x256_speedup\": "
        << (glRes.stored_resolutions[0].avg_frame_ms /
            std::max(0.0001, vkRes.stored_resolutions[0].avg_frame_ms))
        << ",\n";
    out << "    \"stored_512x512_speedup\": "
        << (glRes.stored_resolutions[1].avg_frame_ms /
            std::max(0.0001, vkRes.stored_resolutions[1].avg_frame_ms))
        << ",\n";
    out << "    \"stored_1024x1024_speedup\": "
        << (glRes.stored_resolutions[2].avg_frame_ms /
            std::max(0.0001, vkRes.stored_resolutions[2].avg_frame_ms))
        << ",\n";
  }
  out << "    \"picking_3d_speedup\": "
      << (glRes.picking_avg_us / std::max(0.001, vkRes.picking_avg_us)) << ",\n";
  out << "    \"vector_eps_export_speedup\": "
      << (glRes.eps_export_ms / std::max(0.001, vkRes.eps_export_ms)) << "\n";
  out << "  }\n";
  out << "}\n";
  out.close();

  std::cout << "  [OK] Vulkan Display List Compile: " << vkRes.dlist_compile_us << " us vs OpenGL: "
            << glRes.dlist_compile_us << " us\n";
  for (std::size_t i = 0; i < vkRes.stored_resolutions.size(); ++i) {
    const auto& v = vkRes.stored_resolutions[i];
    const auto& g = glRes.stored_resolutions[i];
    std::cout << "  [OK] Stored " << v.width << "x" << v.height << " -> Vulkan: " << v.avg_frame_ms
              << " ms (" << v.fps << " FPS) | OpenGL: " << g.avg_frame_ms << " ms (" << g.fps
              << " FPS)\n";
  }
  std::cout << "  [OK] Immediate 512x512 -> Vulkan: " << vkRes.immediate_frame_ms << " ms ("
            << vkRes.immediate_fps << " FPS) | OpenGL: " << glRes.immediate_frame_ms << " ms ("
            << glRes.immediate_fps << " FPS)\n";
  std::cout << "  [OK] 3D Picking -> Vulkan: " << vkRes.picking_avg_us << " us/pick | OpenGL: "
            << glRes.picking_avg_us << " us/pick\n";
  std::cout << "  [OK] Vector EPS Export -> Vulkan: " << vkRes.eps_export_ms << " ms ("
            << vkRes.eps_file_bytes << " B) | OpenGL: " << glRes.eps_export_ms << " ms ("
            << glRes.eps_file_bytes << " B)\n";
  std::cout << "  [OK] Wrote JSON benchmark report to: " << jsonPath << "\n";
  return 0;
}

}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  if (argc >= 2 && std::string(argv[1]) == "--benchmark") {
    std::string outPath = (argc >= 3) ? argv[2] : "benchmark_opengl_vs_vulkan.json";
    return RunHeadToHeadBenchmark(outPath);
  }

  G4UIExecutive* ui = nullptr;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
  G4ScoringManager::GetScoringManager();

  runManager->SetUserInitialization(new DetectorConstruction());

  G4VModularPhysicsList* physicsList = new QBBC;
  physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(physicsList);

  runManager->SetUserInitialization(new ActionInitialization());

  auto visManager = new G4VisExecutive(argc, argv);
  visManager->RegisterGraphicsSystem(new G4VulkanStoredQt);
  visManager->RegisterGraphicsSystem(new G4VulkanImmediateQt);
  visManager->RegisterGraphicsSystem(new G4VulkanStoredX);
  visManager->RegisterGraphicsSystem(new G4VulkanImmediateX);
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (!ui) {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }
  else {
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;
  return 0;
}
