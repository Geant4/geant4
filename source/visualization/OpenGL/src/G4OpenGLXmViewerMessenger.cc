//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4OpenGLXmViewerMessenger.cc,v 1.2 2005/11/24 10:23:43 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmViewerMessenger.hh"

#include "G4OpenGLXmViewer.hh"
#include "G4OpenGLXmSliderBar.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"

G4OpenGLXmViewerMessenger::G4OpenGLXmViewerMessenger
(G4OpenGLXmViewer* viewer, const G4String& viewerShortName):
  fpViewer(viewer),
  fViewerShortName(viewerShortName)
{
  G4String directoryName;
  G4String guidance;
  G4String commandName;
  G4bool omitable;

  directoryName = "/vis/oglxm-" + fViewerShortName + '/';
  fpDirectory = new G4UIdirectory(directoryName);
  guidance = "G4OpenGLXmViewer \"" + fViewerShortName + "\" commands.";
  fpDirectory -> SetGuidance(guidance);

  directoryName = directoryName + "set/";
  fpDirectorySet = new G4UIdirectory (directoryName);
  guidance = "G4OpenGLXmViewer \"" + fViewerShortName + "\" set commands.";
  fpDirectorySet -> SetGuidance(guidance);

  commandName = directoryName + "dolly-high";
  fpCommandSetDollyHigh = new G4UIcmdWithADoubleAndUnit
  (commandName, this);
  fpCommandSetDollyHigh -> SetGuidance("Higher limit of dolly slider.");
  fpCommandSetDollyHigh -> SetParameterName("dolly-high", omitable = false);

  commandName = directoryName + "dolly-low";
  fpCommandSetDollyLow = new G4UIcmdWithADoubleAndUnit
  (commandName, this);
  fpCommandSetDollyLow -> SetGuidance("Lower limit of dolly slider.");
  fpCommandSetDollyLow -> SetParameterName("dolly-low", omitable = false);

  commandName = directoryName + "pan-high";
  fpCommandSetPanHigh = new G4UIcmdWithADoubleAndUnit
  (commandName, this);
  fpCommandSetPanHigh -> SetGuidance("Higher limit of pan slider.");
  fpCommandSetPanHigh -> SetParameterName("pan-high", omitable = false);

  commandName = directoryName + "rotation-high";
  fpCommandSetRotationHigh = new G4UIcmdWithADoubleAndUnit
  (commandName, this);
  fpCommandSetRotationHigh -> SetGuidance("Higher limit of rotation slider.");
  fpCommandSetRotationHigh -> SetParameterName("rotation-high", omitable = false);

  commandName = directoryName + "zoom-high";
  fpCommandSetZoomHigh = new G4UIcmdWithADouble
  (commandName, this);
  fpCommandSetZoomHigh -> SetGuidance("Higher limit of zoom slider.");
  fpCommandSetZoomHigh -> SetParameterName("zoom-high", omitable = false);

  commandName = directoryName + "zoom-low";
  fpCommandSetZoomLow = new G4UIcmdWithADouble
  (commandName, this);
  fpCommandSetZoomLow -> SetGuidance("Lower limit of zoom slider.");
  fpCommandSetZoomLow -> SetParameterName("zoom-low", omitable = false);
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
  G4bool panningControlPanel = true;
  G4bool rotationControlPanel = true;

  if (command == fpCommandSetDollyHigh)
    {
      if (fpViewer->fpdolly_slider)
	{
	  fpViewer->dolly_high =
	    fpCommandSetDollyHigh->GetNewDoubleValue(newValue);
	  fpViewer->fpdolly_slider->SetMaxValue (fpViewer->dolly_high);
	  if (fpViewer->fVP.GetDolly() > fpViewer->dolly_high)
	    {
	      fpViewer->fpdolly_slider->SetInitialValue (fpViewer->dolly_high);
	      fpViewer->fVP.SetDolly(fpViewer->dolly_high);
	    }
	  else
	    {
	      fpViewer->fpdolly_slider->SetInitialValue (fpViewer->fVP.GetDolly());
	    }
	}
      else
	{
	  panningControlPanel = false;
	}
    }

  else if (command == fpCommandSetDollyLow)
    {
      if (fpViewer->fpdolly_slider)
	{
	  fpViewer->dolly_low =
	    fpCommandSetDollyLow->GetNewDoubleValue(newValue);
	  fpViewer->fpdolly_slider->SetMinValue (fpViewer->dolly_low);
	  if (fpViewer->fVP.GetDolly() < fpViewer->dolly_low)
	    {
	      fpViewer->fpdolly_slider->SetInitialValue (fpViewer->dolly_low);
	      fpViewer->fVP.SetDolly(fpViewer->dolly_low);
	    }
	  else
	    {
	      fpViewer->fpdolly_slider->SetInitialValue (fpViewer->fVP.GetDolly());
	    }
	}
      else
	{
	  panningControlPanel = false;
	}
    }

  else if (command == fpCommandSetPanHigh)
    {
      if (fpViewer->fppanning_slider)
	{
	  fpViewer->pan_sens_limit =
	    fpCommandSetPanHigh->GetNewDoubleValue(newValue);
	  fpViewer->fppanning_slider->SetMaxValue (fpViewer->pan_sens_limit);
	  fpViewer->fppanning_slider->SetInitialValue (fpViewer->pan_sens_limit / 2.);
	}
      else
	{
	  panningControlPanel = false;
	}
    }

  else if (command == fpCommandSetRotationHigh)
    {
      if (fpViewer->fprotation_slider)
	{
	  // Internally in OpenGLXm, it's in degrees...
	  fpViewer->rot_sens_limit =
	    fpCommandSetRotationHigh->GetNewDoubleValue(newValue) / deg;
	  fpViewer->fprotation_slider->SetMaxValue (fpViewer->rot_sens_limit);
	  fpViewer->fprotation_slider->SetInitialValue (fpViewer->rot_sens_limit / 2.);
	}
      else
	{
	  rotationControlPanel = false;
	}
    }

  else if (command == fpCommandSetZoomHigh)
    {
      if (fpViewer->fpzoom_slider)
	{
	  fpViewer->zoom_high =
	    fpCommandSetZoomHigh->GetNewDoubleValue(newValue);
	  fpViewer->fpzoom_slider->SetMaxValue (fpViewer->zoom_high);
	  fpViewer->fpzoom_slider->SetInitialValue (fpViewer->fVP.GetZoomFactor());
	  if (fpViewer->fVP.GetZoomFactor() > fpViewer->zoom_high)
	    {
	      fpViewer->fpzoom_slider->SetInitialValue (fpViewer->zoom_high);
	      fpViewer->fVP.SetZoomFactor(fpViewer->zoom_high);
	    }
	  else
	    {
	      fpViewer->fpzoom_slider->SetInitialValue (fpViewer->fVP.GetZoomFactor());
	    }
	}
      else
	{
	  panningControlPanel = false;
	}
    }

  else if (command == fpCommandSetZoomLow)
    {
      if (fpViewer->fpzoom_slider)
	{
	  fpViewer->zoom_low =
	    fpCommandSetZoomLow->GetNewDoubleValue(newValue);
	  fpViewer->fpzoom_slider->SetMinValue (fpViewer->zoom_low);
	  fpViewer->fpzoom_slider->SetInitialValue (fpViewer->fVP.GetZoomFactor());
	  if (fpViewer->fVP.GetZoomFactor() < fpViewer->zoom_low)
	    {
	      fpViewer->fpzoom_slider->SetInitialValue (fpViewer->zoom_low);
	      fpViewer->fVP.SetZoomFactor(fpViewer->zoom_low);
	    }
	  else
	    {
	      fpViewer->fpzoom_slider->SetInitialValue (fpViewer->fVP.GetZoomFactor());
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
	    "*** G4OpenGLXmViewerMessenger::SetNewValue: pull down panning"
	    "\n*** control panel and re-issue command."
		 << G4endl;
	  return;
    }

  if (!rotationControlPanel)
    {
	  G4cout << 
	    "*** G4OpenGLXmViewerMessenger::SetNewValue: pull down rotation"
	    "\n*** control panel and re-issue command."
		 << G4endl;
	  return;
    }

  G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/rebuild");
}

#endif
