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
// Ported to Vulkan from the OpenGL work by Andrew Walkden (16th April 1997).
// G4VulkanXmStyleCallbacks :
//                       Several callback functions used by
//                       elements of the control panel to
//                       determine how to visualize the view.

#include "G4VulkanXmViewer.hh"
#include "G4SystemOfUnits.hh"

void G4VulkanXmViewer::drawing_style_callback(Widget w, XtPointer clientData, XtPointer)
{
  G4long choice = (G4long)clientData;
  G4VulkanXmViewer* pView;
  XtVaGetValues(XtParent(w), XmNuserData, &pView, NULL);
  G4ViewParameters::DrawingStyle style;

  switch (choice)
  {
    case 0:
      style = G4ViewParameters::wireframe;
      break;

    case 1:
      style = G4ViewParameters::hlr;
      break;

    case 2:
      style = G4ViewParameters::hsr;
      break;

    case 3:
      style = G4ViewParameters::hlhsr;
      break;

    default:
      style = G4ViewParameters::wireframe;
      G4Exception("G4VulkanXmViewer::drawing_style_callback", "vulkan2006", FatalException,
                  "Unrecognised case in drawing_style_callback.");
  }

  pView->fVP.SetDrawingStyle(style);

  pView->SetView();
  pView->ClearView();
  pView->DrawView();
}

void G4VulkanXmViewer::background_color_callback(Widget w, XtPointer clientData, XtPointer)
{
  G4long choice = (G4long)clientData;
  G4VulkanXmViewer* pView;
  XtVaGetValues(XtParent(w), XmNuserData, &pView, NULL);

  // I need to revisit the kernel if the background colour changes and
  // hidden line removal is enabled, because hlr drawing utilises the
  // background colour in its drawing...
  //  (Note added by JA 13/9/2005) Background now handled in view
  //  parameters.  A kernel visit is triggered on change of background.
  switch (choice)
  {
    case 0:
      ((G4ViewParameters&)pView->GetViewParameters())
        .SetBackgroundColour(G4Colour(1., 1., 1.));  // White
      break;

    case 1:
      ((G4ViewParameters&)pView->GetViewParameters())
        .SetBackgroundColour(G4Colour(0., 0., 0.));  // Black
      break;

    default:
      G4Exception("G4VulkanXmViewer::background_color_callback", "vulkan2008", FatalException,
                  "Unrecognised case in background_color_callback.");
  }

  pView->SetView();
  pView->ClearView();
  pView->DrawView();
}

void G4VulkanXmViewer::transparency_callback(Widget w, XtPointer clientData, XtPointer)
{
  G4long choice = (G4long)clientData;
  G4VulkanXmViewer* pView;
  XtVaGetValues(XtParent(w), XmNuserData, &pView, NULL);

  switch (choice)
  {
    case 0:
      pView->transparency_enabled = false;
      break;

    case 1:
      pView->transparency_enabled = true;
      break;

    default:
      G4Exception("G4VulkanXmViewer::transparency_callback", "vulkan2009", FatalException,
                  "Unrecognised case in transparency_callback.");
  }

  pView->SetNeedKernelVisit(true);
  pView->SetView();
  pView->ClearView();
  pView->DrawView();
}

void G4VulkanXmViewer::antialias_callback(Widget w, XtPointer clientData, XtPointer)
{
  G4long choice = (G4long)clientData;
  G4VulkanXmViewer* pView;
  XtVaGetValues(XtParent(w), XmNuserData, &pView, NULL);

  switch (choice)
  {
    case 0:
      pView->antialiasing_enabled = false;
      vkDisable(VK_LINE_SMOOTH);
      vkDisable(VK_POLYGON_SMOOTH);
      break;

    case 1:
      pView->antialiasing_enabled = true;
      vkEnable(VK_LINE_SMOOTH);
      vkHint(VK_LINE_SMOOTH_HINT, VK_NICEST);
      vkEnable(VK_POLYGON_SMOOTH);
      vkHint(VK_POLYGON_SMOOTH_HINT, VK_NICEST);
      break;

    default:
      G4Exception("G4VulkanXmViewer::antialias_callback", "vulkan2010", FatalException,
                  "Unrecognised case in antialiasing_callback.");
  }

  pView->SetView();
  pView->ClearView();
  pView->DrawView();
}

void G4VulkanXmViewer::haloing_callback(Widget w, XtPointer clientData, XtPointer)
{
  G4long choice = (G4long)clientData;
  G4VulkanXmViewer* pView;
  XtVaGetValues(XtParent(w), XmNuserData, &pView, NULL);

  switch (choice)
  {
    case 0:
      pView->haloing_enabled = false;
      break;

    case 1:
      pView->haloing_enabled = true;
      break;

    default:
      G4Exception("G4VulkanXmViewer::haloing_callback", "vulkan2011", FatalException,
                  "Unrecognised case in haloing_callback.");
  }

  pView->SetView();
  pView->ClearView();
  pView->DrawView();
}

void G4VulkanXmViewer::aux_edge_callback(Widget w, XtPointer clientData, XtPointer)
{
  G4long choice = (G4long)clientData;
  G4VulkanXmViewer* pView;
  XtVaGetValues(XtParent(w), XmNuserData, &pView, NULL);

  switch (choice)
  {
    case 0:
      pView->fVP.SetAuxEdgeVisible(false);
      break;

    case 1:
      pView->fVP.SetAuxEdgeVisible(true);
      break;

    default:
      G4Exception("G4VulkanXmViewer::aux_edge_callback", "vulkan2012", FatalException,
                  "Unrecognised case in aux_edge_callback.");
  }

  pView->SetNeedKernelVisit(true);
  pView->SetView();
  pView->ClearView();
  pView->DrawView();
}

void G4VulkanXmViewer::projection_callback(Widget w, XtPointer clientData, XtPointer)
{
  G4VulkanXmViewer* pView = (G4VulkanXmViewer*)clientData;

  G4int choice = get_int_userData(w);

  switch (choice)
  {
    case 0: {
      pView->fVP.SetFieldHalfAngle(0.);
      break;
    }

    case 1: {
      if (pView->fov > 89.5 || pView->fov <= 0.0)
      {
        G4cout << "Field half angle should be 0 < angle <= 89.5 degrees.";
        G4cout << G4endl;
      }
      else
      {
        pView->fVP.SetFieldHalfAngle(pView->fov * deg);
      }
      break;
    }
    default: {
      G4Exception("G4VulkanXmViewer::projection_callback", "vulkan2013", FatalException,
                  "Unrecognised choice made in projection_callback");
    }
  }

  pView->SetView();
  pView->ClearView();
  pView->DrawView();
}
