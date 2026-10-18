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
// G4VulkanXViewer : Class to provide XWindows specific
//                 functionality for Vulkan in GEANT4

#include "G4VulkanXViewer.hh"

#include "G4LogicalVolume.hh"
#include "G4Normal3D.hh"
#include "G4VulkanFontBaseStore.hh"
#include "G4VulkanSceneHandler.hh"
#include "G4Point3D.hh"
#include "G4StateManager.hh"
#include "G4Text.hh"
#include "G4Threading.hh"
#include "G4VSolid.hh"
#include "G4VisExtent.hh"
#include "G4VisManager.hh"

#include <X11/Xatom.h>
#include <X11/Xmu/StdCmap.h>
#include <X11/Xutil.h>
#include <assert.h>

#include <chrono>
#include <sstream>
#include <thread>

int G4VulkanXViewer::snglBuf_RGBA[12] = {1, 8, 1, 8, 1, 8, 1, 24, 1, 8, 1, None};

int G4VulkanXViewer::dblBuf_RGBA[13] = {1, 8, 1, 8, 1, 8, 1, 1, 24, 1, 8, 1, None};

#define NewString(str)                                                                     \
  ((str) != 0                                                                              \
     ? (strncpy((char*)malloc((unsigned)strlen(str) + 1), str, (unsigned)strlen(str) + 1)) \
     : (char*)0)

#define USE_DEFAULT_COLORMAP 1
#define USE_STANDARD_COLORMAP 0

XVisualInfo* G4VulkanXViewer::vi_single_buffer = 0;
XVisualInfo* G4VulkanXViewer::vi_double_buffer = 0;

extern "C"
{
  static Bool G4VulkanXViewerWaitForNotify(Display*, XEvent* e, char* arg)
  {
    return (e->type == MapNotify) && (e->xmap.window == (Window)arg);
  }
}

void G4VulkanXViewer::PresentVulkanFramebuffer()
{
  if (!dpy || !win || !vi) return;
  auto& ctx = G4VulkanContext::Instance();
  int w = ctx.fFbWidth;
  int h = ctx.fFbHeight;
  if (w <= 0 || h <= 0 || ctx.fColorBuffer.size() < static_cast<size_t>(w * h * 4)) return;
  char* imgData = static_cast<char*>(std::malloc(static_cast<size_t>(w) * h * 4));
  if (!imgData) return;
  for (int y = 0; y < h; ++y)
  {
    const uint8_t* srcRow = ctx.fColorBuffer.data() + static_cast<size_t>(h - 1 - y) * w * 4;
    uint8_t* dstRow = reinterpret_cast<uint8_t*>(imgData) + static_cast<size_t>(y) * w * 4;
    for (int x = 0; x < w; ++x)
    {
      dstRow[x * 4 + 0] = srcRow[x * 4 + 2];
      dstRow[x * 4 + 1] = srcRow[x * 4 + 1];
      dstRow[x * 4 + 2] = srcRow[x * 4 + 0];
      dstRow[x * 4 + 3] = 255;
    }
  }
  XImage* ximg = XCreateImage(dpy, vi->visual, vi->depth, ZPixmap, 0, imgData, w, h, 32, w * 4);
  if (ximg)
  {
    GC gc = XCreateGC(dpy, win, 0, nullptr);
    XPutImage(dpy, win, gc, ximg, 0, 0, 0, 0, w, h);
    XFreeGC(dpy, gc);
    XDestroyImage(ximg);
  }
  else
  {
    std::free(imgData);
  }
  XFlush(dpy);
}

void G4VulkanXViewer::SetView()
{
  G4VulkanContext::Instance().EnsureVulkanInitialised();
  G4VulkanViewer::SetView();
}

void G4VulkanXViewer::ShowView()
{
  vkFlush();
  PresentVulkanFramebuffer();

  if (fVP.IsPicking())
  {
    G4cout << "Window activated for picking (left-mouse), exit (middle-mouse)." << G4endl;
    while (true)
    {
      if (XPending(dpy))
      {
        XNextEvent(dpy, &event);
        if (event.type == ButtonPress && event.xbutton.button == 1)
        {
          G4cout << Pick(event.xbutton.x, event.xbutton.y) << G4endl;
        }
        else if (event.type == ButtonPress && event.xbutton.button == 2)
          break;
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
  }
}

void G4VulkanXViewer::SwitchToVisSubThread()
{
#ifdef G4MULTITHREADED
  G4VulkanContext::Instance().EnsureVulkanInitialised();
  cxVisSubThread = G4VulkanContext::Instance().fInstance;
#endif
}

void G4VulkanXViewer::SwitchToMasterThread()
{
#ifdef G4MULTITHREADED
  G4VulkanContext::Instance().EnsureVulkanInitialised();
  cxVisSubThread = VK_NULL_HANDLE;
#endif
}

void G4VulkanXViewer::GetXConnection()
{
  // get a connection.
  dpy = XOpenDisplay(0);  // Uses DISPLAY environment variable.
  if (!dpy)
  {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4VulkanXViewer::G4VulkanXViewer couldn't open display." << G4endl;
    return;
  }

  // Ensure Vulkan context is initialized.
  G4VulkanContext::Instance().EnsureVulkanInitialised();
}

void G4VulkanXViewer::CreateVKXContext(XVisualInfo* v)
{
  vi = v;
  // get window's attributes
  if (!XGetWindowAttributes(dpy, XRootWindow(dpy, vi->screen), &xwa))
  {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4VulkanXViewer::G4VulkanXViewer couldn't return window attributes" << G4endl;
    return;
  }

  // initialize the master Vulkan instance handle
  G4VulkanContext::Instance().EnsureVulkanInitialised();
  cxMaster = G4VulkanContext::Instance().fInstance;

  // New stab at getting a colormap

  Status status;
  int i, numCmaps;

  status = XmuLookupStandardColormap(dpy, vi->screen, vi->visualid, vi->depth, XA_RGB_BEST_MAP,
                                     False, True);

  if (status == 1)
  {
    cmap = 0;
    XStandardColormap* standardCmaps = XAllocStandardColormap();
    status = XGetRGBColormaps(dpy, XRootWindow(dpy, vi->screen), &standardCmaps, &numCmaps,
                              XA_RGB_BEST_MAP);
    if (status == 1)
    {
      for (i = 0; i < numCmaps; i++)
      {
        if (standardCmaps[i].visualid == vi->visualid)
        {
          cmap = standardCmaps[i].colormap;
          break;
        }
      }
    }
    XFree(standardCmaps);
    if (cmap)
    {
      if (G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
        G4cout << "Got standard cmap" << G4endl;
    }
    else
    {
      // if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      //	G4cerr << "G4VulkanXViewer::G4VulkanXViewer failed to allocate a standard colormap."
      //	       << G4endl;
      cmap = XCreateColormap(dpy, XRootWindow(dpy, vi->screen), vi->visual, AllocNone);
      if (cmap)
      {
        if (G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
          G4cout << "Created own cmap" << G4endl;
      }
      // G.Barrand : at end, we should do a XFreeColormap(dpy,cmap) when cmap is no more used.
    }
  }
  else
  {
    cmap = XCreateColormap(dpy, XRootWindow(dpy, vi->screen), vi->visual, AllocNone);
    if (cmap)
    {
      if (G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
        G4cout << "Created own cmap" << G4endl;
    }
    // G.Barrand : at end, we should do a XFreeColormap(dpy,cmap) when cmap is no more used.
  }

  if (!cmap)
  {
    fViewId = -1;  // This flags an error.
    if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
      G4cout << "G4VulkanXViewer::G4VulkanXViewer failed to allocate a Colormap." << G4endl;
    return;
  }
}

void G4VulkanXViewer::CreateMainWindow()
{
  // create a window
  swa.colormap = cmap;
  swa.border_pixel = 0;
  swa.event_mask = ExposureMask | ButtonPressMask | StructureNotifyMask;
  swa.backing_store = WhenMapped;

  // Window size and position...
  size_hints = XAllocSizeHints();

  ResizeWindow(fVP.GetWindowSizeHintX(), fVP.GetWindowSizeHintY());

  G4int x_origin = fVP.GetWindowAbsoluteLocationHintX(DisplayWidth(dpy, vi->screen));

  // FIXME,  screen size != window size on MAC, but I don't know have to get the menuBar
  // size on MAC. L.Garnier 01/2009
  G4int y_origin = fVP.GetWindowAbsoluteLocationHintY(DisplayHeight(dpy, vi->screen));

  size_hints->base_width = getWinWidth();
  size_hints->base_height = getWinHeight();
  size_hints->x = x_origin;
  size_hints->y = y_origin;
  if (fVP.IsWindowSizeHintX() && fVP.IsWindowLocationHintX() && fVP.IsWindowLocationHintY())
  {
    size_hints->flags |= PSize | PPosition;
  }
  else if (fVP.IsWindowSizeHintX() && !(fVP.IsWindowLocationHintX() || fVP.IsWindowLocationHintY()))
  {
    size_hints->flags |= PSize;
  }
  else if ((!fVP.IsWindowSizeHintX()) && fVP.IsWindowLocationHintX() && fVP.IsWindowLocationHintY())
  {
    size_hints->flags |= PPosition;
  }
  if (G4VisManager::GetVerbosity() >= G4VisManager::confirmations)
    G4cout << "Window name: " << fName << G4endl;
  strncpy(charViewName, fName, 99);
  charViewName[99] = '\0';
  char* window_name = charViewName;
  char* icon_name = charViewName;
  // char tmpatom[] = "XA_WM_NORMAL_HINTS";
  wm_hints = XAllocWMHints();
  class_hints = XAllocClassHint();

  XStringListToTextProperty(&window_name, 1, &windowName);
  XStringListToTextProperty(&icon_name, 1, &iconName);

  wm_hints->initial_state = NormalState;
  wm_hints->input = True;
  wm_hints->icon_pixmap = icon_pixmap;
  wm_hints->flags = StateHint | IconPixmapHint | InputHint;

  class_hints->res_name = NewString("G4Vulkan");
  class_hints->res_class = NewString("G4Vulkan");

  win = XCreateWindow(dpy, XRootWindow(dpy, vi->screen), x_origin, y_origin, getWinWidth(),
                      getWinHeight(), 0, vi->depth, InputOutput, vi->visual,
                      CWBorderPixel | CWColormap | CWEventMask | CWBackingStore, &swa);

  XSetWMProperties(dpy, win, &windowName, &iconName, 0, 0, size_hints, wm_hints, class_hints);

  // request X to Draw window on screen.
  XMapWindow(dpy, win);

  // Wait for window to appear (wait for an "expose" event).
  XIfEvent(dpy, &event, G4VulkanXViewerWaitForNotify, (char*)win);

  // Initialize Vulkan context for window
  G4VulkanContext::Instance().EnsureVulkanInitialised();
}

void G4VulkanXViewer::CreateFontLists()
{
  std::map<G4double, G4String> fonts;  // G4VMarker screen size and font name.
  fonts[10.] = "-adobe-courier-bold-r-normal--10-100-75-75-m-60-iso8859-1";
  fonts[11.] = "-adobe-courier-bold-r-normal--11-80-100-100-m-60-iso8859-1";
  fonts[12.] = "-adobe-courier-bold-r-normal--12-120-75-75-m-70-iso8859-1";
  fonts[13.] = "fixed";
  fonts[14.] = "-adobe-courier-bold-r-normal--14-100-100-100-m-90-iso8859-1";
  fonts[17.] = "-adobe-courier-bold-r-normal--17-120-100-100-m-100-iso8859-1";
  fonts[18.] = "-adobe-courier-bold-r-normal--18-180-75-75-m-110-iso8859-1";
  fonts[20.] = "-adobe-courier-bold-r-normal--20-140-100-100-m-110-iso8859-1";
  fonts[24.] = "-adobe-courier-bold-r-normal--24-240-75-75-m-150-iso8859-1";
  fonts[25.] = "-adobe-courier-bold-r-normal--25-180-100-100-m-150-iso8859-1";
  fonts[34.] = "-adobe-courier-bold-r-normal--34-240-100-100-m-200-iso8859-1";
  std::map<G4double, G4String>::const_iterator i;
  for (i = fonts.begin(); i != fonts.end(); ++i)
  {
    XFontStruct* font_info = XLoadQueryFont(dpy, i->second);
    if (!font_info)
    {
      G4cerr << "G4VulkanXViewer::CreateFontLists XLoadQueryFont failed for font\n  " << i->second
             << G4endl;
      continue;
    }
    G4int font_base = vkGenLists(256);
    if (!font_base)
    {
      G4cerr << "G4VulkanXViewer::CreateFontLists out of display lists for fonts." << G4endl;
      XFreeFont(dpy, font_info);
      continue;
    }
    G4int width = font_info->max_bounds.width;
    G4VulkanFontBaseStore::AddFontBase(this, font_base, i->first, i->second, width);
    XFreeFont(dpy, font_info);
  }
}

void G4VulkanXViewer::DrawText(const G4Text& g4text)
{
  if (isVk2psWriting())
  {
    G4VulkanViewer::DrawText(g4text);
  }
  else
  {
    G4VSceneHandler::MarkerSizeType sizeType;
    G4double size = fSceneHandler.GetMarkerSize(g4text, sizeType);

    const G4VulkanFontBaseStore::FontInfo& fontInfo =
      G4VulkanFontBaseStore::GetFontInfo(this, (int)size);
    if (fontInfo.fFontBase < 0)
    {
      static G4int callCount = 0;
      ++callCount;
      // if (callCount <= 10 || callCount%100 == 0) {
      if (callCount <= 1)
      {
        G4cout << "G4VulkanXViewer::DrawText: No fonts available for \"" << fName
               << "\"\n  Called with " << g4text << G4endl;
      }
      return;
    }

    const G4Colour& c = fSceneHandler.GetTextColour(g4text);
    vkColor4d(c.GetRed(), c.GetGreen(), c.GetBlue(), c.GetAlpha());

    G4Point3D position = g4text.GetPosition();

    G4String textString = g4text.GetText();
    const char* textCString = textString.c_str();

    // Set position for raster-style drawers (X, Xm)
    vkRasterPos3d(position.x(), position.y(), position.z());

    vkPushAttrib(VK_LIST_BIT);

    // Calculate move for centre and right adjustment
    G4double span = textString.size() * fontInfo.fWidth;
    G4double xmove = 0., ymove = 0.;
    switch (g4text.GetLayout())
    {
      case G4Text::left:
        break;
      case G4Text::centre:
        xmove -= span / 2.;
        break;
      case G4Text::right:
        xmove -= span;
    }

    // Add offsets
    xmove += g4text.GetXOffset();
    ymove += g4text.GetYOffset();

    // Do move
    vkBitmap(0, 0, 0, 0, xmove, ymove, 0);

    // Write characters
    vkListBase(fontInfo.fFontBase);
    vkCallLists((G4int)strlen(textCString), VK_UNSIGNED_BYTE, (VKubyte*)textCString);
    vkPopAttrib();
  }
}

G4VulkanXViewer::G4VulkanXViewer(G4VulkanSceneHandler& scene)
  : G4VViewer(scene, -1), G4VulkanViewer(scene), vi_immediate(0), vi_stored(0), vi(0), cmap(0)
{
  // To satisfy Coverity
  xwa.visual = 0;
  iconName.value = 0;
  xwa.screen = 0;
  windowName.value = 0;

  GetXConnection();
  if (fViewId < 0) return;

  // Select a TrueColor X11 visual for Vulkan presentation
  if (!vi_single_buffer)
  {
    XVisualInfo vTemplate;
    vTemplate.screen = XDefaultScreen(dpy);
    vTemplate.depth = DefaultDepth(dpy, XDefaultScreen(dpy));
    vTemplate.c_class = TrueColor;
    int numVisuals = 0;
    vi_single_buffer = XGetVisualInfo(dpy, VisualScreenMask | VisualDepthMask | VisualClassMask,
                                      &vTemplate, &numVisuals);
    if (!vi_single_buffer)
    {
      vi_single_buffer = XGetVisualInfo(dpy, VisualScreenMask, &vTemplate, &numVisuals);
    }
  }
  if (!vi_double_buffer)
  {
    vi_double_buffer = vi_single_buffer;
  }

  if (vi_single_buffer || vi_double_buffer)
  {
    if (!vi_double_buffer)
    {
      G4cout << "G4VulkanXViewer::G4VulkanXViewer: unable to get a double buffer visual."
                "\n  Working with a single buffer."
             << G4endl;
    }
  }
  else
  {
    if (!vi_single_buffer)
    {
      G4cout << "G4VulkanXViewer::G4VulkanXViewer: unable to get a single buffer visual." << G4endl;
    }
    if (!vi_double_buffer)
    {
      G4cout << "G4VulkanXViewer::G4VulkanXViewer: unable to get a double buffer visual." << G4endl;
    }
  }

  if (vi_single_buffer)
  {
    vi_immediate = vi_single_buffer;
    attributeList = snglBuf_RGBA;
  }

  if (!vi_immediate)
  {
    // next try for a double buffered RGB, but Draw to top buffer
    if (vi_double_buffer)
    {
      vi_immediate = vi_double_buffer;
      attributeList = dblBuf_RGBA;
    }
  }

  // Now try for a visual suitable for VulkanStored...
  // Try for a double buffered RGB window
  if (vi_double_buffer)
  {
    vi_stored = vi_double_buffer;
    attributeList = dblBuf_RGBA;
  }

  if (!vi_immediate || !vi_stored)
  {
    G4cout << "G4VulkanXViewer::G4VulkanXViewer: unable to get required visuals." << G4endl;
    fViewId = -1;  // This flags an error.
  }

  //  vkClearColor (0., 0., 0., 0.);
  //  vkClearDepth (1.);
}

G4VulkanXViewer::~G4VulkanXViewer()
{
  if (fViewId >= 0)
  {
    // Close a window from here
    cxMaster = VK_NULL_HANDLE;
    if (win) XDestroyWindow(dpy, win);  // ...if already deleted in
    // sub-class G4VulkanXmViewer.
    // We should do a XFreeColormap(dpy,cmap); if cmap had been get with XCreateColormap.
    XFlush(dpy);
  }
}
