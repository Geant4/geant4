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

#include "G4VtkOffscreenViewer.hh"
#include "G4VtkSceneHandler.hh"
#include "G4VSceneHandler.hh"

#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"

G4VtkOffscreenViewer::G4VtkOffscreenViewer(G4VSceneHandler& sceneHandler, const G4String& name)
  : G4VtkViewer(sceneHandler, name)
{
}

G4VtkOffscreenViewer::~G4VtkOffscreenViewer()
{
}

void G4VtkOffscreenViewer::Initialise()
{
  G4VtkViewer::Initialise();
  _renderWindow->SetOffScreenRendering(1);
}

void G4VtkOffscreenViewer::FinishView()
{
  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  fVtkSceneHandler.Modified();

  _renderWindow->Render();

  vtkNew<vtkWindowToImageFilter> windowToImageFilter;
  windowToImageFilter->SetInput(_renderWindow);
  windowToImageFilter->Update();

  vtkNew<vtkPNGWriter> writer;
  writer->SetFileName("offscreen.png");
  writer->SetInputConnection(windowToImageFilter->GetOutputPort());
  writer->Write();
}
