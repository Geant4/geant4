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
// $Id: G4OpenGLXmViewerMessenger.cc 66373 2012-12-18 09:41:34Z gcosmo $

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmViewerMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4OpenGLXmViewer.hh"
#include "G4OpenGLXmSliderBar.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"

#include "G4VisManager.hh"

G4OpenGLXmViewerMessenger* G4OpenGLXmViewerMessenger::fpInstance = 0;

G4OpenGLXmViewerMessenger* G4OpenGLXmViewerMessenger::GetInstance()
{
  if (!fpInstance) fpInstance = new G4OpenGLXmViewerMessenger;
  return fpInstance;
}

G4OpenGLXmViewerMessenger::G4OpenGLXmViewerMessenger()
{
  G4bool omitable;

  fpDirectory = new G4UIdirectory("/vis/oglxm/");
  fpDirectory->SetGuidance("G4OpenGLXmViewer commands.");

  fpDirectorySet = new G4UIdirectory ("/vis/oglxm/set/");
  fpDirectorySet->SetGuidance("G4OpenGLXmViewer set commands.");

  fpCommandSetDollyHigh =
    new G4UIcmdWithADoubleAndUnit("/vis/oglxm/set/dolly-high", this);
  fpCommandSetDollyHigh->SetGuidance("Higher limit of dolly slider.");
  fpCommandSetDollyHigh->SetParameterName("dolly-high", omitable = false);

  fpCommandSetDollyLow =
    new G4UIcmdWithADoubleAndUnit("/vis/oglxm/set/dolly-low", this);
  fpCommandSetDollyLow->SetGuidance("Lower limit of dolly slider.");
  fpCommandSetDollyLow->SetParameterName("dolly-low", omitable = false);

  fpCommandSetPanHigh =
    new G4UIcmdWithADoubleAndUnit("/vis/oglxm/set/pan-high", this);
  fpCommandSetPanHigh->SetGuidance("Higher limit of pan slider.");
  fpCommandSetPanHigh->SetParameterName("pan-high", omitable = false);

  fpCommandSetRotationHigh =
    new G4UIcmdWithADoubleAndUnit("/vis/oglxm/set/rotation-high", this);
  fpCommandSetRotationHigh->SetGuidance("Higher limit of rotation slider.");
  fpCommandSetRotationHigh->SetParameterName("rotation-high", omitable = false);

  fpCommandSetZoomHigh =
    new G4UIcmdWithADouble("/vis/oglxm/set/zoom-high", this);
  fpCommandSetZoomHigh->SetGuidance("Higher limit of zoom slider.");
  fpCommandSetZoomHigh->SetParameterName("zoom-high", omitable = false);

  fpCommandSetZoomLow =
    new G4UIcmdWithADouble("/vis/oglxm/set/zoom-low", this);
  fpCommandSetZoomLow->SetGuidance("Lower limit of zoom slider.");
  fpCommandSetZoomLow->SetParameterName("zoom-low", omitable = false);
}

G4OpenGLXmViewerMessenger::~G4OpenGLXmViewerMessenger ()
{
  delete fpCommandSetZoomLow;
  delete fpCommandSetZoomHigh;
  delete fpCommandSetRotationHigh;
  delete fpCommandSetPanHigh;
  delete fpCommandSetDollyLow;
  delete fpCommandSetDollyHigh;
  delete fpDirectorySet;
  delete fpDirectory;
}

void G4OpenGLXmViewerMessenger::SetNewValue
(G4UIcommand* command, G4String newValue)
{
  G4VisManager* pVisManager = G4VisManager::GetInstance();

  G4VViewer* pVViewer = pVisManager->GetCurrentViewer();

  if (!pVViewer) {
    G4cout <<
      "G4OpenGLXmViewerMessenger::SetNewValue: No current viewer."
      "\n  \"/vis/open\", or similar, to get one."
           << G4endl;
    return;
  }

  G4OpenGLXmViewer* pViewer = dynamic_cast<G4OpenGLXmViewer*>(pVViewer);

  if (!pViewer) {
    G4cout <<
      "G4OpenGLXmViewerMessenger::SetNewValue: Current viewer is not of type"
      "\n  OGLIXm or OGLSXm.  Use \"/vis/viewer/select\" or \"/vis/open\"."
           << G4endl;
    return;
  }

  G4bool panningControlPanel = true;
  G4bool rotationControlPanel = true;

  if (command == fpCommandSetDollyHigh)
    {
      if (pViewer->fpdolly_slider)
	{
	  pViewer->dolly_high =
	    fpCommandSetDollyHigh->GetNewDoubleValue(newValue);
	  pViewer->fpdolly_slider->SetMaxValue (pViewer->dolly_high);
	  if (pViewer->fVP.GetDolly() > pViewer->dolly_high)
	    {
	      pViewer->fpdolly_slider->SetInitialValue (pViewer->dolly_high);
	      pViewer->fVP.SetDolly(pViewer->dolly_high);
	    }
	  else
	    {
	      pViewer->fpdolly_slider->SetInitialValue (pViewer->fVP.GetDolly());
	    }
	}
      else
	{
	  panningControlPanel = false;
	}
    }

  else if (command == fpCommandSetDollyLow)
    {
      if (pViewer->fpdolly_slider)
	{
	  pViewer->dolly_low =
	    fpCommandSetDollyLow->GetNewDoubleValue(newValue);
	  pViewer->fpdolly_slider->SetMinValue (pViewer->dolly_low);
	  if (pViewer->fVP.GetDolly() < pViewer->dolly_low)
	    {
	      pViewer->fpdolly_slider->SetInitialValue (pViewer->dolly_low);
	      pViewer->fVP.SetDolly(pViewer->dolly_low);
	    }
	  else
	    {
	      pViewer->fpdolly_slider->SetInitialValue (pViewer->fVP.GetDolly());
	    }
	}
      else
	{
	  panningControlPanel = false;
	}
    }

  else if (command == fpCommandSetPanHigh)
    {
      if (pViewer->fppanning_slider)
	{
	  pViewer->pan_sens_limit =
	    fpCommandSetPanHigh->GetNewDoubleValue(newValue);
	  pViewer->fppanning_slider->SetMaxValue (pViewer->pan_sens_limit);
	  pViewer->fppanning_slider->SetInitialValue (pViewer->pan_sens_limit / 2.);
	}
      else
	{
	  panningControlPanel = false;
	}
    }

  else if (command == fpCommandSetRotationHigh)
    {
      if (pViewer->fprotation_slider)
	{
	  // Internally in OpenGLXm, it's in degrees...
	  pViewer->rot_sens_limit =
	    fpCommandSetRotationHigh->GetNewDoubleValue(newValue) / deg;
	  pViewer->fprotation_slider->SetMaxValue (pViewer->rot_sens_limit);
	  pViewer->fprotation_slider->SetInitialValue (pViewer->rot_sens_limit / 2.);
	}
      else
	{
	  rotationControlPanel = false;
	}
    }

  else if (command == fpCommandSetZoomHigh)
    {
      if (pViewer->fpzoom_slider)
	{
	  pViewer->zoom_high =
	    fpCommandSetZoomHigh->GetNewDoubleValue(newValue);
	  pViewer->fpzoom_slider->SetMaxValue (pViewer->zoom_high);
	  pViewer->fpzoom_slider->SetInitialValue (pViewer->fVP.GetZoomFactor());
	  if (pViewer->fVP.GetZoomFactor() > pViewer->zoom_high)
	    {
	      pViewer->fpzoom_slider->SetInitialValue (pViewer->zoom_high);
	      pViewer->fVP.SetZoomFactor(pViewer->zoom_high);
	    }
	  else
	    {
	      pViewer->fpzoom_slider->SetInitialValue (pViewer->fVP.GetZoomFactor());
	    }
	}
      else
	{
	  panningControlPanel = false;
	}
    }

  else if (command == fpCommandSetZoomLow)
    {
      if (pViewer->fpzoom_slider)
	{
	  pViewer->zoom_low =
	    fpCommandSetZoomLow->GetNewDoubleValue(newValue);
	  pViewer->fpzoom_slider->SetMinValue (pViewer->zoom_low);
	  pViewer->fpzoom_slider->SetInitialValue (pViewer->fVP.GetZoomFactor());
	  if (pViewer->fVP.GetZoomFactor() < pViewer->zoom_low)
	    {
	      pViewer->fpzoom_slider->SetInitialValue (pViewer->zoom_low);
	      pViewer->fVP.SetZoomFactor(pViewer->zoom_low);
	    }
	  else
	    {
	      pViewer->fpzoom_slider->SetInitialValue (pViewer->fVP.GetZoomFactor());
	    }
	}
      else
	{
	  panningControlPanel = false;
	}
    }

  if (!panningControlPanel)
    {
	  G4cout << 
	    "G4OpenGLXmViewerMessenger::SetNewValue: pull down panning"
	    "\n  control panel and re-issue command."
		 << G4endl;
	  return;
    }

  if (!rotationControlPanel)
    {
	  G4cout << 
	    "G4OpenGLXmViewerMessenger::SetNewValue: pull down rotation"
	    "\n  control panel and re-issue command."
		 << G4endl;
	  return;
    }

  G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/rebuild");
}

#endif
