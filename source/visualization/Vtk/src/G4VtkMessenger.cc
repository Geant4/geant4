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

#include "G4Tokenizer.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4VisManager.hh"
#include "G4VtkSceneHandler.hh"
#include "G4VtkViewer.hh"

#include "vtkObject.h"

G4VtkMessenger* G4VtkMessenger::fpInstance = nullptr;

G4VtkMessenger* G4VtkMessenger::GetInstance()
{
  if (fpInstance == nullptr) fpInstance = new G4VtkMessenger;
  return fpInstance;
}

G4VtkMessenger::G4VtkMessenger()
{
  G4bool omitable;

  /***************************** Vtk directory *****************************/
  fpDirectory = new G4UIdirectory("/vis/vtk/", true);
  fpDirectory->SetGuidance("G4VtkViewer commands.");

  // Clear non-G4
  fpCommandClearNonG4 = new G4UIcommand("/vis/vtk/clearNonG4", this);
  fpCommandClearNonG4->SetGuidance("Clear non G4 objects from visualisation (image/3d overlays");

  // Export command
  fpCommandExport = new G4UIcommand("/vis/vtk/export", this);
  fpCommandExport->SetGuidance("Export a screenshot or OBJ file of current Vtk viewer");

  // File type for export
  auto parameterExport = new G4UIparameter("format", 's', omitable = true);
  fpCommandExport->SetGuidance("File type (jpg,tiff,eps,ps,obj,vtp,gltf,vrml)");
  fpCommandExport->SetParameter(parameterExport);

  // File name for export
  parameterExport = new G4UIparameter("file-name", 's', omitable = true);
  fpCommandExport->SetGuidance("File name");
  fpCommandExport->SetParameter(parameterExport);

  // Export cutter command
  fpCommandExportCutter = new G4UIcommand("/vis/vtk/exportCutter", this);
  fpCommandExportCutter->SetGuidance("Export a VTP file of the cutters if they exist");

  // File name for export
  auto parameterExportCutter = new G4UIparameter("file-name", 's', omitable = true);
  parameterExportCutter->SetGuidance("File name");
  parameterExportCutter->SetDefaultValue("cutter.vtp");
  fpCommandExportCutter->SetParameter(parameterExportCutter);

  // Vtk debug print
  fpCommandDebugPrint = new G4UIcommand("/vis/vtk/printDebug", this);
  fpCommandDebugPrint->SetGuidance("Debug print of pipelines");

  /***************************** Set directory *****************************/
  fpDirectorySet = new G4UIdirectory("/vis/vtk/set/", true);
  fpDirectorySet->SetGuidance("G4VtkViewer set commands");

  // Vtk warnings output
  fpCommandWarnings = new G4UIcmdWithABool("/vis/vtk/set/warnings", this);
  fpCommandWarnings->SetParameterName("enable-warnings", false);
  fpCommandWarnings->SetGuidance("Enable (True) or disable (False) VTK warnings");

  // HUD command
  fpCommandHUD = new G4UIcmdWithABool("/vis/vtk/set/hud", this);
  fpCommandHUD->SetGuidance("Enable or disable HUD for VTK");

  // Clipper command
  fpCommandClipper = new G4UIcommand("/vis/vtk/set/clipper", this);
  fpCommandClipper->SetGuidance("Enable a cutaway plane (clipper) in VTK");
  auto fpCommandClipperParam = new G4UIparameter("widget", 'b', omitable = true);
  fpCommandClipperParam->SetDefaultValue(1);
  fpCommandClipper->SetParameter(fpCommandClipperParam);

  // Clutter command
  fpCommandCutter = new G4UIcommand("/vis/vtk/set/cutter", this);
  fpCommandCutter->SetGuidance("Enable a section plane (cutter) in VTK");
  auto fpCommandCutterParam = new G4UIparameter("widget", 'b', omitable = true);
  fpCommandCutterParam->SetDefaultValue(1);
  fpCommandCutter->SetParameter(fpCommandCutterParam);

  // Vtk polyhedron pipeline selection
  auto fpCommandPolyhedronPipelineParam = new G4UIparameter("format", 's', omitable = true);
  fpCommandPolyhedronPipeline = new G4UIcommand("/vis/vtk/set/polyhedronPipeline", this);
  fpCommandPolyhedronPipeline->SetGuidance("Select type of polyhedron pipeline");
  fpCommandPolyhedronPipeline->SetGuidance("Type (separate, tensor, append, bake)");
  fpCommandPolyhedronPipeline->SetParameter(fpCommandPolyhedronPipelineParam);

  // Shadows command
  fpCommandShadow = new G4UIcommand("/vis/vtk/set/shadows", this);
  fpCommandClipper->SetGuidance("Enable/disable shadows in VTK");
  auto fpCommandShadowParam = new G4UIparameter("enable", 'b', omitable = true);
  fpCommandShadowParam->SetDefaultValue(1);
  fpCommandShadow->SetParameter(fpCommandShadowParam);

  /***************************** Add directory *****************************/
  fpDirectoryAdd = new G4UIdirectory("/vis/vtk/add/", true);
  fpDirectoryAdd->SetGuidance("G4VtkViewer add commands");

  fpCommandImageOverlay = new G4UIcommand("/vis/vtk/add/imageOverlay", this);
  auto parameterImageOverlay = new G4UIparameter("filename", 's', omitable = false);
  parameterImageOverlay->SetGuidance("Image file name ");
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);

  parameterImageOverlay = new G4UIparameter("imageBottomLeftX", 'd', omitable = false);
  parameterImageOverlay->SetGuidance("image bottom left x");
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);

  parameterImageOverlay = new G4UIparameter("imageBottomLeftY", 'd', omitable = false);
  parameterImageOverlay->SetGuidance("image bottom left y");
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);

  parameterImageOverlay = new G4UIparameter("imageTopRightX", 'd', omitable = false);
  parameterImageOverlay->SetGuidance("image top right x");
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);

  parameterImageOverlay = new G4UIparameter("imageTopRightY", 'd', omitable = false);
  parameterImageOverlay->SetGuidance("image top right y");
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);

  parameterImageOverlay = new G4UIparameter("worldBottomLeftX", 'd', omitable = false);
  parameterImageOverlay->SetGuidance("world bottom left x");
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);

  parameterImageOverlay = new G4UIparameter("worldBottomLeftY", 'd', omitable = false);
  parameterImageOverlay->SetGuidance("world bottom left y");
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);

  parameterImageOverlay = new G4UIparameter("worldTopRightX", 'd', omitable = false);
  parameterImageOverlay->SetGuidance("world top right x");
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);

  parameterImageOverlay = new G4UIparameter("worldTopRightY", 'd', omitable = false);
  parameterImageOverlay->SetGuidance("world top right y");
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);

  parameterImageOverlay = new G4UIparameter("xRotation", 'd', omitable = true);
  parameterImageOverlay->SetGuidance("rotation in x");
  parameterImageOverlay->SetDefaultValue(0.0);
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);

  parameterImageOverlay = new G4UIparameter("yRotation", 'd', omitable = true);
  parameterImageOverlay->SetGuidance("rotation in y");
  parameterImageOverlay->SetDefaultValue(0.0);
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);

  parameterImageOverlay = new G4UIparameter("zRotation", 'd', omitable = true);
  parameterImageOverlay->SetGuidance("rotation in z");
  parameterImageOverlay->SetDefaultValue(0.0);
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);

  parameterImageOverlay = new G4UIparameter("xTranslation", 'd', omitable = true);
  parameterImageOverlay->SetGuidance("translation in x");
  parameterImageOverlay->SetDefaultValue(0.0);
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);

  parameterImageOverlay = new G4UIparameter("yTranslation", 'd', omitable = true);
  parameterImageOverlay->SetGuidance("translation in y");
  parameterImageOverlay->SetDefaultValue(0.0);
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);

  parameterImageOverlay = new G4UIparameter("zTranslation", 'd', omitable = true);
  parameterImageOverlay->SetGuidance("translation in z");
  parameterImageOverlay->SetDefaultValue(0.0);
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);

  parameterImageOverlay = new G4UIparameter("alpha", 'd', omitable = true);
  parameterImageOverlay->SetGuidance("alpha");
  parameterImageOverlay->SetDefaultValue(0.5);
  fpCommandImageOverlay->SetParameter(parameterImageOverlay);
}

G4VtkMessenger::~G4VtkMessenger()
{
  delete fpDirectory;
  delete fpDirectorySet;
  delete fpDirectoryAdd;
  delete fpCommandExport;
  delete fpCommandExportCutter;
  delete fpCommandWarnings;
  delete fpCommandDebugPrint;
  delete fpCommandPolyhedronPipeline;
  delete fpCommandImageOverlay;
}

G4String G4VtkMessenger::GetCurrentValue(G4UIcommand* /*command*/)
{
  return G4String();
}

void G4VtkMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  G4VisManager* pVisManager = G4VisManager::GetInstance();

  G4VViewer* pViewer = pVisManager->GetCurrentViewer();
  if (pViewer == nullptr) {
    G4cout << "G4VtkMessenger::SetNewValue: No current viewer.\n"
           << "\"/vis/open\", or similar, to get one." << G4endl;
    return;
  }

  auto* pVtkViewer = dynamic_cast<G4VtkViewer*>(pViewer);
  if (pVtkViewer == nullptr) {
    G4cout << "G4VtkMessenger::SetNewValue: Current viewer is not of type VTK. \n"
           << "(It is \"" << pViewer->GetName() << "\".)\n"
           << "Use \"/vis/viewer/select\" or \"/vis/open\"." << G4endl;
    return;
  }

  if (command == fpCommandClearNonG4) {
    auto sceneHandler = dynamic_cast<G4VtkSceneHandler*>(pViewer->GetSceneHandler());
    auto transientStore = sceneHandler->GetTransientStore();

    transientStore.ClearNonG4();
  }
  else if (command == fpCommandExport) {
    G4String format, name;

    std::istringstream iss(newValue);
    iss >> format >> name;

    if (format == "jpg" || format == "tiff" || format == "png" || format == "bmp" || format == "pnm"
        || format == "ps")
      pVtkViewer->ExportScreenShot(name, format);
    else if (format == "obj")
      pVtkViewer->ExportOBJScene(name);
    else if (format == "vrml")
      pVtkViewer->ExportVRMLScene(name);
    else if (format == "vtp")
      pVtkViewer->ExportVTPScene(name);
    else if (format == "gltf")
      pVtkViewer->ExportGLTFScene(name);
    else
      G4cout << "Unknown /vis/vtk/export file format" << G4endl;
  }
  else if (command == fpCommandExportCutter) {
    std::istringstream iss(newValue);

    G4String fileName;
    iss >> fileName;
    pVtkViewer->ExportVTPCutter(fileName);
  }
  else if (command == fpCommandWarnings) {
    if (G4UIcommand::ConvertToBool(newValue)) {
      vtkObject::GlobalWarningDisplayOn();
    }
    else {
      vtkObject::GlobalWarningDisplayOff();
    }
  }
  else if (command == fpCommandHUD) {
    if (G4UIcommand::ConvertToBool(newValue)) {
      pVtkViewer->EnableHUD();
    }
    else {
      pVtkViewer->DisableHUD();
    }
  }
  else if (command == fpCommandDebugPrint) {
    pVtkViewer->Print();
  }
  else if (command == fpCommandPolyhedronPipeline) {
    G4String temp;

    std::istringstream iss(newValue);

    G4String pipelineType;
    iss >> pipelineType;

    pVtkViewer->SetPolyhedronPipeline(pipelineType);
  }
  else if (command == fpCommandImageOverlay) {
    G4String temp;

    G4String fileName;
    G4double imageBottomLeft[2] = {0, 0};
    G4double imageTopRight[2] = {1, 1};
    G4double worldBottomLeft[2] = {0, 0};
    G4double worldTopRight[2] = {1, 1};
    G4double rotation[3] = {0, 0, 0};
    G4double translation[3] = {0, 0, 0};
    G4double alpha = 0;

    std::istringstream iss(newValue);

    iss >> fileName >> imageBottomLeft[0] >> imageBottomLeft[1] >> imageTopRight[0]
      >> imageTopRight[1] >> worldBottomLeft[0] >> worldBottomLeft[1] >> worldTopRight[0]
      >> worldTopRight[1] >> rotation[0] >> rotation[1] >> rotation[2] >> translation[0]
      >> translation[1] >> translation[2] >> alpha;

    pVtkViewer->AddImageOverlay(fileName, alpha, imageBottomLeft, worldBottomLeft, imageTopRight,
                                worldTopRight, rotation, translation);
  }
  else if (command == fpCommandClipper) {
    pVtkViewer->EnableClipper(G4Plane3D(), true);
  }
  else if (command == fpCommandCutter) {
    pVtkViewer->EnableCutter(G4Plane3D(), true);
  }
  else if (command == fpCommandShadow) {
    pVtkViewer->EnableShadows();
  }
}
