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

#include "G4VtkMessenger.hh"
#include "G4VtkViewer.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"

#include "G4VisManager.hh"

#include "vtkObject.h"

G4VtkMessenger* G4VtkMessenger::fpInstance = nullptr;

G4VtkMessenger* G4VtkMessenger::GetInstance()
{
  if (!fpInstance) fpInstance = new G4VtkMessenger;
  return fpInstance;
}

G4VtkMessenger::G4VtkMessenger() {
  G4bool omitable;

  fpDirectory = new G4UIdirectory("/vis/vtk/");
  fpDirectory->SetGuidance("G4VtkViewer commands.");

  // Export command
  fpCommandExport = new G4UIcommand("/vis/vtk/export", this);
  fpCommandExport->SetGuidance ("Export a screenshot or OBJ file of current Vtk viewer");

  // File type for export
  auto parameterExport = new G4UIparameter ("name", 's', omitable = true);
  fpCommandExport->SetGuidance ("File type (jpg,tiff,eps,ps,obj,vrml)");
  fpCommandExport->SetParameter(parameterExport);

  // File name for export
  parameterExport = new G4UIparameter ("name", 's', omitable = true);
  fpCommandExport->SetGuidance ("File name");
  fpCommandExport->SetParameter(parameterExport);

  // Vtk warnings output
  fpCommandWarnings = new G4UIcmdWithABool("/vis/vtk/warnings", this);
  fpCommandExport->SetGuidance ("Enable (True) or disable (False) VTK warnings");
}

G4VtkMessenger::~G4VtkMessenger()
{
  delete fpDirectory;
  delete fpCommandExport;
  delete fpCommandWarnings;
}
G4String G4VtkMessenger::GetCurrentValue(G4UIcommand * /*command*/) {
  return G4String();
}

void G4VtkMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{

  G4VisManager* pVisManager = G4VisManager::GetInstance();

  G4VViewer* pViewer = pVisManager->GetCurrentViewer();
  if (!pViewer) {
    G4cout << "G4VtkMessenger::SetNewValue: No current viewer.\n"
           << "\"/vis/open\", or similar, to get one."
           << G4endl;
    return;
  }

  auto* pVtkViewer = dynamic_cast<G4VtkViewer*>(pViewer);
  if (!pVtkViewer) {
    G4cout << "G4VtkMessenger::SetNewValue: Current viewer is not of type VTK. \n"
           << "(It is \""
           << pViewer->GetName()
           << "\".)\n"
           << "Use \"/vis/viewer/select\" or \"/vis/open\"."
           << G4endl;
    return;
  }

  if (command == fpCommandExport)
  {
    G4String format, name;

    std::istringstream iss(newValue);
    iss >> format >> name;

    if(format == "jpg" || format == "tiff" ||
       format == "png" || format == "bmp" ||
       format == "pnm" || format == "ps")
      pVtkViewer->ExportScreenShot(name, format);
    else if(format == "obj")
      pVtkViewer->ExportOBJScene(name);
    else if(format == "vrml")
      pVtkViewer->ExportVRMLScene(name);
    else if(format == "vtp")
      pVtkViewer->ExportVTPScene(name);
    else
      G4cout << "Unknown /vis/vtk/export file format" << G4endl;
  }
  else if (command == fpCommandWarnings)
  {
    vtkObject::GlobalWarningDisplayOff();
  }
}
