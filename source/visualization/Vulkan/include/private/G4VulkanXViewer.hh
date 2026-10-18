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
//                   functionality for Vulkan in GEANT4

#ifndef G4VULKANXVIEWER_HH
#define G4VULKANXVIEWER_HH

#include "G4Vulkan.hh"
#include "G4VulkanViewer.hh"
#include "globals.hh"

#include <X11/Intrinsic.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

class G4VulkanSceneHandler;
class G4Text;

class G4VulkanXViewer : virtual public G4VulkanViewer
{
    friend class G4VulkanXViewerMessenger;
    friend class G4VulkanXmViewer;

  public:

    G4VulkanXViewer(G4VulkanSceneHandler& scene);
    virtual ~G4VulkanXViewer();
    void SetView();
    void ShowView();
    void SwitchToVisSubThread();
    void SwitchToMasterThread();
    void DrawText(const G4Text&);

  protected:

    void GetXConnection();
    void CreateVKXContext(XVisualInfo* vi);
    void PresentVulkanFramebuffer();
    virtual void CreateMainWindow();
    virtual void CreateFontLists();

    static int snglBuf_RGBA[12];
    static int dblBuf_RGBA[13];

    //////////////////////////////Pixmap (screen dump) production functions/////

    XWindowAttributes xwa;
    Display* dpy;
    static XVisualInfo* vi_single_buffer;
    static XVisualInfo* vi_double_buffer;
    XVisualInfo *vi_immediate, *vi_stored, *vi;
    Colormap cmap;
    XSetWindowAttributes swa;
    Window win;
    VkInstance cxMaster;
    VkInstance cxVisSubThread;
    XEvent event;
    G4int *attributeList, errorBase, eventBase, major, minor;
    XSizeHints* norm_hints;
    XWMHints* wm_hints;
    XClassHint* class_hints;
    Pixmap icon_pixmap;
    XSizeHints* size_hints;
    Atom Xatom;
    XTextProperty windowName, iconName;
    char charViewName[100];

  private:

    G4VulkanXViewer(const G4VulkanXViewer&);
    G4VulkanXViewer& operator=(const G4VulkanXViewer&);
    VkInstance tmp_cx;
};

#endif
